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
! CONTRIBUTIONS: Marco Bacci                                               !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
!---------------------------------------------------------------------
!
! ROUTINES NEEDED FOR SOLVENT-ACCESSIBLE VOLUME CALCULATIONS
!
!---------------------------------------------------------------------
!
!
subroutine setup_freesolv()
!
  use iounit
  use aminos
  use sequen
  use fos
  use fyoc
!
  implicit none
!
  character(3) resname
  integer i,lk
!
! note that this is our first sanity check whether a residue is actually
! supported. while parameters might be read in, this gives a change to
! bail out, if implementation isn't there yet
!
  do i=1,nseq
    resname = amino(seqtyp(i))
    if (disulf(i).gt.0) then
      lk = crlk_idx(i)
      if ((crosslink(lk)%itstype.eq.1).OR.(crosslink(lk)%itstype.eq.2)) then
!       disulfide linkage
        freesolv(i,1:3) = FOS_LINK_CC(1:3)
      else
        write(ilog,*) 'Fatal. Encountered unsupported crosslink type in setup_freesolv(...).&
 & This is most likely an omission bug.'
        call fexit()
      end if
      cycle
    end if
!
!   residues that just have backbone are skipped (they're still
!   supported of course), unsupported ones require patch
    if ((resname.eq.'GLY').OR.(resname.eq.'ACE').OR.&
 &      (resname.eq.'NME').OR.(resname.eq.'FOR').OR.&
 &      (resname.eq.'NH2').OR.(resname.eq.'DIB').OR.&
 &      (resname.eq.'D5P').OR.(resname.eq.'RIB').OR.&
 &      (resname.eq.'R5P').OR.(resname.eq.'UNK')) then
      freesolv(i,1:3) = 0.0
    else if (resname.eq.'ALA') then
      freesolv(i,1:3) = FOS_ALA(1:3)
    else if (resname.eq.'AIB') then
      freesolv(i,1:3) = FOS_AIB(1:3)
    else if (resname.eq.'ABA') then
      freesolv(i,1:3) = FOS_ABA(1:3)
    else if (resname.eq.'VAL') then
      freesolv(i,1:3) = FOS_VAL(1:3)
    else if (resname.eq.'NVA') then
      freesolv(i,1:3) = FOS_NVA(1:3)
    else if (resname.eq.'LEU') then
      freesolv(i,1:3) = FOS_LEU(1:3)
    else if (resname.eq.'ILE') then
      freesolv(i,1:3) = FOS_ILE(1:3)
    else if (resname.eq.'NLE') then
      freesolv(i,1:3) = FOS_NLE(1:3)
    else if (resname.eq.'MET') then
      freesolv(i,1:3) = FOS_MET(1:3)
    else if (resname.eq.'SER') then
      freesolv(i,1:3) = FOS_SER(1:3)
    else if (resname.eq.'SEP') then
      freesolv(i,1:3) = FOS_SEP(1:3)
    else if (resname.eq.'THR') then
      freesolv(i,1:3) = FOS_THR(1:3)
    else if (resname.eq.'TPO') then
      freesolv(i,1:3) = FOS_TPO(1:3)
    else if (resname.eq.'CYS') then
      freesolv(i,1:3) = FOS_CYS(1:3)
    else if (resname.eq.'CYX') then
      freesolv(i,1:3) = FOS_CYX(1:3)
    else if (resname.eq.'PHE') then
      freesolv(i,1:3) = FOS_PHE(1:3)
    else if (resname.eq.'TYR') then
      freesolv(i,1:3) = FOS_TYR(1:3)
    else if (resname.eq.'TYO') then
      freesolv(i,1:3) = FOS_TYO(1:3)
    else if (resname.eq.'PTR') then
      freesolv(i,1:3) = FOS_PTR(1:3)
    else if (resname.eq.'TRP') then
      freesolv(i,1:3) = FOS_TRP(1:3)
    else if (resname.eq.'ASN') then
      freesolv(i,1:3) = FOS_ASN(1:3)
    else if (resname.eq.'GLN') then
      freesolv(i,1:3) = FOS_GLN(1:3)
    else if (resname.eq.'ASP') then
      freesolv(i,1:3) = FOS_ASP(1:3)
    else if (resname.eq.'ASH') then
      freesolv(i,1:3) = FOS_ASH(1:3)
    else if (resname.eq.'GLU') then
      freesolv(i,1:3) = FOS_GLU(1:3)
    else if (resname.eq.'GLH') then
      freesolv(i,1:3) = FOS_GLH(1:3)
    else if (resname.eq.'ARG') then
      freesolv(i,1:3) = FOS_ARG(1:3)
    else if (resname.eq.'LYS') then
      freesolv(i,1:3) = FOS_LYS(1:3)
    else if (resname.eq.'LYD') then
      freesolv(i,1:3) = FOS_LYD(1:3)
    else if (resname.eq.'KAC') then
      freesolv(i,1:3) = FOS_KAC(1:3)
    else if (resname.eq.'KM1') then
      freesolv(i,1:3) = FOS_KM1(1:3)
    else if (resname.eq.'KM2') then
      freesolv(i,1:3) = FOS_KM2(1:3)
    else if (resname.eq.'KM3') then
      freesolv(i,1:3) = FOS_KM3(1:3)
    else if (resname.eq.'PRO') then
      freesolv(i,1:3) = FOS_PRO(1:3)
    else if (resname.eq.'HYP') then
      freesolv(i,1:3) = FOS_HYP(1:3)
    else if (resname.eq.'HID') then
      freesolv(i,1:3) = FOS_HID(1:3)
    else if (resname.eq.'HIE') then
      freesolv(i,1:3) = FOS_HIE(1:3)
    else if (resname.eq.'HIP') then
      freesolv(i,1:3) = FOS_HIP(1:3)
    else if (resname.eq.'ORN') then
      freesolv(i,1:3) = FOS_ORN(1:3)
    else if (resname.eq.'DAB') then
      freesolv(i,1:3) = FOS_DAB(1:3)
    else if (resname.eq.'O2 ') then
      freesolv(i,1:3) = FOS_O2(1:3)
    else if (resname.eq.'SPC') then
      freesolv(i,1:3) = FOS_SPC(1:3)
    else if (resname.eq.'T3P') then
      freesolv(i,1:3) = FOS_T3P(1:3)
    else if (resname.eq.'NA+') then
      freesolv(i,1:3) = FOS_NA(1:3)
    else if (resname.eq.'CL-') then
      freesolv(i,1:3) = FOS_CL(1:3)
    else if (resname.eq.'K+ ') then
      freesolv(i,1:3) = FOS_K(1:3)
    else if (resname.eq.'BR-') then
      freesolv(i,1:3) = FOS_BR(1:3)
    else if (resname.eq.'CS+') then
      freesolv(i,1:3) = FOS_CS(1:3)
    else if (resname.eq.'I- ') then
      freesolv(i,1:3) = FOS_I(1:3)
    else if (resname.eq.'NH4') then
      freesolv(i,1:3) = FOS_NH4(1:3)
    else if (resname.eq.'GDN') then
      freesolv(i,1:3) = FOS_GDN(1:3)
    else if (resname.eq.'AC-') then
      freesolv(i,1:3) = FOS_AC(1:3)
    else if (resname.eq.'1MN') then
      freesolv(i,1:3) = FOS_1MN(1:3)
    else if (resname.eq.'2MN') then
      freesolv(i,1:3) = FOS_2MN(1:3)
    else if (resname.eq.'LCP') then
      freesolv(i,1:3) = FOS_LCP(1:3)
    else if (resname.eq.'NO3') then
      freesolv(i,1:3) = FOS_NO3(1:3)
    else if (resname.eq.'URE') then
      freesolv(i,1:3) = FOS_URE(1:3)
    else if (resname.eq.'NMF') then
      freesolv(i,1:3) = FOS_NMF(1:3)
    else if (resname.eq.'FOA') then
      freesolv(i,1:3) = FOS_FOA(1:3)
    else if (resname.eq.'CH4') then
      freesolv(i,1:3) = FOS_CH4(1:3)
    else if (resname.eq.'PRP') then
      freesolv(i,1:3) = FOS_PRP(1:3)
    else if (resname.eq.'IBU') then
      freesolv(i,1:3) = FOS_IBU(1:3)
    else if (resname.eq.'NBU') then
      freesolv(i,1:3) = FOS_NBU(1:3)
    else if (resname.eq.'MSH') then
      freesolv(i,1:3) = FOS_MSH(1:3)
    else if (resname.eq.'EOH') then
      freesolv(i,1:3) = FOS_EOH(1:3)
    else if (resname.eq.'EMT') then
      freesolv(i,1:3) = FOS_EMT(1:3)
    else if (resname.eq.'BEN') then
      freesolv(i,1:3) = FOS_BEN(1:3)
    else if (resname.eq.'NAP') then
      freesolv(i,1:3) = FOS_NAP(1:3)
    else if (resname.eq.'TOL') then
      freesolv(i,1:3) = FOS_TOL(1:3)
    else if (resname.eq.'IMD') then
      freesolv(i,1:3) = FOS_IMD(1:3)
    else if (resname.eq.'IME') then
      freesolv(i,1:3) = FOS_IME(1:3)
    else if (resname.eq.'MIN') then
      freesolv(i,1:3) = FOS_MIN(1:3)
    else if (resname.eq.'NMA') then
      freesolv(i,1:3) = FOS_NMA(1:3)
    else if (resname.eq.'ACA') then
      freesolv(i,1:3) = FOS_ACA(1:3)
    else if (resname.eq.'PPA') then
      freesolv(i,1:3) = FOS_PPA(1:3)
    else if (resname.eq.'MOH') then
      freesolv(i,1:3) = FOS_MOH(1:3)
    else if (resname.eq.'PCR') then
      freesolv(i,1:3) = FOS_PCR(1:3)
    else if (resname.eq.'DMA') then
      freesolv(i,1:3) = FOS_DMA(1:3)
    else if (resname.eq.'THY') then
      freesolv(i,1:3) = FOS_THY(1:3)
    else if (resname.eq.'PUR') then
      freesolv(i,1:3) = FOS_PUR(1:3)
    else if (resname.eq.'ADE') then
      freesolv(i,1:3) = FOS_ADE(1:3)
    else if (resname.eq.'GUA') then
      freesolv(i,1:3) = FOS_GUA(1:3)
    else if (resname.eq.'CYT') then
      freesolv(i,1:3) = FOS_CYT(1:3)
    else if (resname.eq.'URA') then
      freesolv(i,1:3) = FOS_URA(1:3)
    else if ((resname.eq.'DPT').OR.(resname.eq.'RPT').OR.&
 &           (resname.eq.'DIT').OR.(resname.eq.'RIT')) then
      freesolv(i,1:3) = FOS_XPT(1:3)
    else if ((resname.eq.'DPA').OR.(resname.eq.'RPA').OR.&
 &           (resname.eq.'DIA').OR.(resname.eq.'RIA')) then
      freesolv(i,1:3) = FOS_XPA(1:3)
    else if ((resname.eq.'DPG').OR.(resname.eq.'RPG').OR.&
 &           (resname.eq.'DIG').OR.(resname.eq.'RIG')) then
      freesolv(i,1:3) = FOS_XPG(1:3)
    else if ((resname.eq.'DPC').OR.(resname.eq.'RPC').OR.&
 &           (resname.eq.'DIC').OR.(resname.eq.'RIC')) then
      freesolv(i,1:3) = FOS_XPC(1:3)
    else if ((resname.eq.'DPU').OR.(resname.eq.'RPU').OR.&
 &           (resname.eq.'DIU').OR.(resname.eq.'RIU')) then
      freesolv(i,1:3) = FOS_XPU(1:3)
    else
      write(ilog,*) 'Fatal. Implicit solvent does not support residu&
 &e ',resname,' at the moment. Check back later.'
      call fexit()
    end if
  end do
!
end
!
!---------------------------------------------------------------------------
!
subroutine get_sexv()
!
  use iounit
!
  implicit none
!
  RTYPE sexv
!
  sexv = 0.0
  write(ilog,*) 'Fatal. Called get_sexv() which is currently an unsupported routine.&
 & This is most certainly an omission bug.'
  call fexit()
!
! unsupported dummy fxn
!
end
!
!---------------------------------------------------------------------
!
subroutine init_svte(isvmode)
!
  use atoms
  use iounit
  use sequen
  use polypep
  use molecule
!
  implicit none
!
  integer i,isvmode,rs,ii
!
  if (isvmode.eq.0) then
    do i=1,n
      svte(i) = 0.0
      svbu(i) = 0.0
    end do
  else if (isvmode.eq.1) then
    do i=1,n
      svbu(i) = svte(i)
      svte(i) = 0.0
    end do
  else if (isvmode.eq.2) then
    do rs=1,nseq
      do i=1,at(rs)%nbb+at(rs)%nsc
        if (i.le.at(rs)%nbb) then
          ii = at(rs)%bb(i)
        else
           ii = at(rs)%sc(i-at(rs)%nbb)
        end if
        svte(ii) = svbu(ii) - svte(ii)
        if (svte(ii).ne.0.0) then
          rs_vec(rs) = 1
          if (which_fosg(ii).gt.100) rs_vec2(rs+1) = 1
        end if
      end do
    end do
  else if (isvmode.eq.3) then
    do i=1,n
      svte(i) = atstv(i)
    end do
  else if (isvmode.eq.4) then
    do i=1,n
      atsav(i) = atsav(i) - svte(i)
    end do
  else if (isvmode.eq.5) then
    do i=1,n
      atsav(i) = atsav(i) + svte(i)
    end do
  else if (isvmode.eq.6) then
    do i=1,n
      atsav(i) = svte(i)
    end do
  else
    write(ilog,*) 'Fatal. Called init_svte(...) with bad mode (',&
 &isvmode,').'
    call fexit()
  end if
!
end
!
!-------------------------------------------------------------------------------
!
#ifdef ENABLE_THREADS
!
subroutine init_svte_threads(isvmode,tpi)
!
  use atoms
  use iounit
  use sequen
  use threads
  use cutoffs, ONLY: rsinfo
!
  implicit none
!
  integer, INTENT(IN):: tpi,isvmode
!
  integer rs,i,sta,sto,stx,sta2,sto2
  RTYPE saveps
!
  saveps = 1.0e-9
  if (tpi.gt.0) then
    saveps = thrdat%maxn*saveps
    if (isvmode.eq.2) then
      sta2 = thr_limits(1,tpi)
      sto2 = thr_limits(2,tpi)
      if (thr_limits(4,tpi).ge.thr_limits(3,tpi)) then
        sta = rsinfo(thr_limits(3,tpi),1)
        sto = rsinfo(thr_limits(4,tpi),1) +  rsinfo(thr_limits(4,tpi),2)
      else
        sta = 1
        sto = 0
      end if
      stx = thr_limits(4,tpi)
    else
      sta = thr_limits(1,tpi)
      sto = thr_limits(2,tpi)
      stx = thr_limits(4,tpi)
    end if
  else
    sta = 1
    sto = n
    stx = nseq + 1
  end if
!
  if (isvmode.eq.0) then
    svte(sta:sto) = 0.0
    svbu(sta:sto) = 0.0
    if (tpi.gt.0) thr_svte(1:n,tpi) = 0.0
  else if (isvmode.eq.1) then
    if (tpi.gt.0) then
      do i=1,thrdat%maxn
        svbu(sta:sto) = svbu(sta:sto) + thr_svte(sta:sto,i)
      end do
!$OMP BARRIER
      thr_svte(1:n,tpi) = 0.0
    else
      svbu(1:n) = svte(1:n)
      svte(1:n) = 0.0
    end if
  else if (isvmode.eq.2) then
    if (tpi.gt.0) then
      do i=1,thrdat%maxn
        svte(sta2:sto2) = svte(sta2:sto2) + thr_svte(sta2:sto2,i)
      end do
    end if
!$OMP BARRIER
    do i=sta,sto
      svte(i) = svbu(i) - svte(i)
      if (abs(svte(i)).gt.saveps) then
        rs = atmres(i)
        rs_vec(rs) = 1
        if (which_fosg(i).gt.100) then
          if (rs.eq.stx) then
!$OMP CRITICAL(SVTE_RSV2)
            rs_vec2(rs+1) = 1
!$OMP END CRITICAL(SVTE_RSV2)
          else
            rs_vec2(rs+1) = 1
          end if
        end if
      end if
    end do
!$OMP BARRIER
  else if (isvmode.eq.3) then
    svte(sta:sto) = atstv(sta:sto)
  else if (isvmode.eq.4) then
    atsav(sta:sto) = atsav(sta:sto) - svte(sta:sto)
  else if (isvmode.eq.5) then
    atsav(sta:sto) = atsav(sta:sto) + svte(sta:sto)
  else if (isvmode.eq.6) then
    if (tpi.gt.0) then
      do i=1,thrdat%maxn
        svte(sta:sto) = svte(sta:sto) + thr_svte(sta:sto,i)
      end do
    end if
!$OMP BARRIER
    atsav(sta:sto) = svte(sta:sto)
  else
    write(ilog,*) 'Fatal. Called init_svte_threads(...) with bad mode (',isvmode,').'
    call fexit()
  end if
!
end
!
#endif
!
!-----------------------------------------------------------------------------------------
!
function ln_ipol(ati)
!
  use atoms
  use energies
!
  implicit none
!
  RTYPE ln_ipol,rat
  integer ati
!
  rat = atsav(ati)/atbvol(ati)
  if (rat.gt.atsavmaxfr(ati)) then
    ln_ipol = 1.0
  else if (rat.gt.par_IMPSOLV(5)) then
    ln_ipol = (1.0/(atsavmaxfr(ati)-par_IMPSOLV(5)))*&
 &                      (rat - par_IMPSOLV(5))
  else
    ln_ipol = 0.0
  end if
!  if (ati.eq.1) then
!  if((ati.eq.3).OR.(ati.eq.11).OR.(ati.eq.10).OR.(ati.eq.17)) then
!    write(ilog,*) ati,' ',ln_ipol
!  end if
!
  return
!
end
!
!
function ln_ipol2(ati)
!
  use atoms
  use energies
!
  implicit none
!
  RTYPE ln_ipol2,rat
  integer ati
!
  rat = (atsav(ati)+svte(ati))/atbvol(ati)
  if (rat.gt.atsavmaxfr(ati)) then
    ln_ipol2 = 1.0
  else if (rat.gt.par_IMPSOLV(5)) then
    ln_ipol2 = (1.0/(atsavmaxfr(ati)-par_IMPSOLV(5)))*&
 &                      (rat - par_IMPSOLV(5))
  else
    ln_ipol2 = 0.0
  end if
!  if (ati.eq.1) then
!  if((ati.eq.3).OR.(ati.eq.11).OR.(ati.eq.10).OR.(ati.eq.17)) then
!    write(ilog,*) ati,' ',ln_ipol2
!  end if
!
  return
!
end
!
!------------------------------------------------------------------------------------------------
!
subroutine do_sav()
!
  use mcsums
  use energies
  use atoms
  use mpistuff
  use molecule
  use sequen
  use system
  use grandensembles
  use pdb
  use fos
!
  implicit none
!
  RTYPE tsv,fs_ipol,dummy,incr
  integer i,imol,b1,b2
  logical sayno,ismember
!
  sayno = .false.
!
  if (use_frame_weights.EQV..true.) then
    incr = 1.0*framewts(curframe)
  else
    incr = 1.0
  end if
!
! if we do not use the potential relying on local densities, we have to compute them first
! the (fixed) topology part is always taken care of beforehand
  if (use_IMPSOLV.EQV..false.) then
    call setup_savol(sayno,tsv)
  end if
!
  tsv = 0.0
  call get_molsav(tsv,dummy,sayno)
!
! average values for total and for individual atoms
  nsavavg = nsavavg + 1
  if (use_frame_weights.EQV..false.) then
    savavg = ((nsavavg-1.0)*savavg + tsv)/(1.0*nsavavg)
  end if
  do imol=1,nmol
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
        if (ismember(ispresent,imol).EQV..false.) cycle
      end if
    end if
    natsavavg(imol) = natsavavg(imol) + 1
    do i=atmol(imol,1),atmol(imol,2)
      atsavavg(1,i) = ((natsavavg(imol)-1.0)*atsavavg(1,i) + incr*atsav(i))/&
 &                  (1.0*natsavavg(imol))
      atsavavg(2,i) = ((natsavavg(imol)-1.0)*atsavavg(2,i) + incr*fs_ipol(i))/&
 &                  (1.0*natsavavg(imol))
    end do
  end do
  if (savreq%nats.gt.0) then
    do i=1,savreq%nats
      b1 = max(1,min(100,ceiling(atsav(savreq%idx(i))/(atbvol(savreq%idx(i))*0.01))))
      savreq%hists(i,1,b1) = savreq%hists(i,1,b1) + incr
      b2 = max(1,min(100,ceiling(fs_ipol(savreq%idx(i))/0.01)))
      savreq%hists(i,2,b2) = savreq%hists(i,2,b2) + incr
    end do
  end if
!
 47   format(i10,10000000(1x,g14.7))
!
  if (mod(nsavavg,savreq%instfreq).eq.0) then
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then
      if (savreq%nats.gt.0) then
        write(isasa,47) nstep,savavg,tsv,atsav(savreq%idx(1:savreq%nats))/atbvol(savreq%idx(1:savreq%nats))
      else
        write(isasa,47) nstep,savavg,tsv
      end if
    else
      if (myrank.eq.0) then
        call MPI_AVGWriteSAV()
      else
        call MPI_AVGSendSAV()
      end if
      if (mpi_cnt_sav.eq.(mpi_nodes-1)) then
        mpi_cnt_sav = 0
      else
        mpi_cnt_sav = mpi_cnt_sav + 1
      end if
    end if
#else
    if (savreq%nats.gt.0) then
      if (use_frame_weights.EQV..false.) then
        write(isasa,47) nstep,savavg,tsv,atsav(savreq%idx(1:savreq%nats))/atbvol(savreq%idx(1:savreq%nats))
      else
        write(isasa,47) nstep,tsv,atsav(savreq%idx(1:savreq%nats))/atbvol(savreq%idx(1:savreq%nats))
      end if
    else
      if (use_frame_weights.EQV..false.) then
        write(isasa,47) nstep,savavg,tsv
      else
        write(isasa,47) nstep,tsv
      end if
    end if
#endif
  end if
!
end
!
!---------------------------------------------------------------------
!
subroutine prt_sav()
!
  use atoms
  use energies
  use mpistuff
  use system
  use molecule
  use pdb
  use iounit
  use fos
  use sequen
!
  implicit none
!
  integer i,k,freeunit,iu,ii,jj,imol
  character(60) fn
#ifdef ENABLE_MPI
  character(re_aux(10)) nod
#endif
  logical exists
  RTYPE, ALLOCATABLE:: normer(:)
  RTYPE tmz
!
 47   format(i10,1x,g14.7,1x,g14.7,1x,g14.7,1x,g14.7,1x,g14.7,1x,g14.7)
 48   format(f8.3,1x,10000000(g13.6,1x,g13.6,1x))
 49   format('#',8x,10000000(10x,i8,10x))
!
  allocate(normer(nmol))
  if (use_frame_weights.EQV..true.) then
    normer(:) = 0.0
    k = 0
    do i=1,framecnt
      if ((framelst2(i).gt.nequil).AND.(mod(framelst2(i),savcalc).eq.0)) then
        k = k + 1
        normer(:) = normer(:) + framewts2(i)
      end if
    end do
    if (k.ne.natsavavg(1)) then
      write(ilog,*) 'Warning. Solvent-accessible volume analyses have inconsistent normalization. This is &
 &a bug. Please report this error and ignore all output in SAV_BY_ATOM.dat, and related files.'
    end if
  else
    normer(:) = 1.0*natsavavg(:)
  end if
!
  if (use_frame_weights.EQV..true.) then
    do imol=1,nmol
      atsavavg(1,atmol(imol,1):atmol(imol,2)) = natsavavg(imol)*atsavavg(1,atmol(imol,1):atmol(imol,2))/normer(imol)
      atsavavg(2,atmol(imol,1):atmol(imol,2)) = natsavavg(imol)*atsavavg(2,atmol(imol,1):atmol(imol,2))/normer(imol)
    end do
  end if
!
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    call int2str(myrank,nod,re_aux(10))
    fn = 'N_'//nod//'_SAV_BY_ATOM.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'SAV_BY_ATOM.dat'
  end if
#else
  fn = 'SAV_BY_ATOM.dat'
#endif
  call strlims(fn,ii,jj)
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    iu = freeunit()
    open(unit=iu,file=fn(ii:jj),status='old')
    close(unit=iu,status='delete')
  end if
  iu=freeunit()
  open(unit=iu,file=fn(ii:jj),status='new')
  do i=1,n
    write(iu,47) i,atsavavg(1,i)/atbvol(i),atsavavg(2,i),&
 &atsavmaxfr(i),atsavprm(i,2),atsavprm(i,5),atsavred(i)
  end do
  close(unit=iu)
!
! finally the histograms for specific atoms
  if (savreq%nats.gt.0) then
!
    exists = .false.
    do k=1,savreq%nats
      if (natsavavg(molofrs(atmres(savreq%idx(k)))).gt.0) then
        exists = .true.
        exit
      end if
    end do
    if (exists.EQV..false.) return
!
    do k=1,savreq%nats
      do i=1,2
        tmz = sum(savreq%hists(k,i,1:100))
        if (tmz.gt.0.0) then
          savreq%hists(k,i,1:100) = savreq%hists(k,i,1:100)/tmz
        end if
      end do
    end do
!
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then
      call int2str(myrank,nod,re_aux(10))
      fn = 'N_'//nod//'_SAV_HISTS.dat'
    else if (use_MPIAVG.EQV..true.) then
      fn = 'SAV_HISTS.dat'
    end if
#else
    fn = 'SAV_HISTS.dat'
#endif
!
    call strlims(fn,ii,jj)
    inquire(file=fn(ii:jj),exist=exists)
    if(exists) then
      iu = freeunit()
      open(unit=iu,file=fn(ii:jj),status='old')
      close(unit=iu,status='delete')
    end if
    iu=freeunit()
    open(unit=iu,file=fn(ii:jj),status='new')
    write(iu,49) savreq%idx(1:savreq%nats)
    do i=1,100
      write(iu,48) (i-0.5)*0.01,(savreq%hists(k,1,i),savreq%hists(k,2,i),k=1,savreq%nats)
    end do
    close(unit=iu)
  end if
! 
end
!
!----------------------------------------------------------------------------------------------------
!
subroutine get_molsav(tsv,tsv2,do2)
!
  use atoms
  use energies
  use math
  use aminos
  use polypep
  use sequen
  use molecule
  use grandensembles
  use system
!
  implicit none
!
  integer imol,rs,ii,i
  character(3) resname
  logical do2,ismember
  RTYPE tsv,tsv2,ln_ipol
  RTYPE svr,svr2,sphrad,cylrad,cylvol,sphvol,ln_ipol2
!
!
  do imol=1,nmol
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
        if (ismember(ispresent,imol).EQV..false.) cycle
      end if
    end if
!   accumulate approximate heavy atom residue-SAV ...
    do rs=rsmol(imol,1),rsmol(imol,2)
      svr = 0.0
      svr2 = 0.0
      resname = amino(seqtyp(rs))
      do i=1,at(rs)%nbb+at(rs)%nsc
        if (i.le.at(rs)%nbb) then
          ii = at(rs)%bb(i)
        else
          ii = at(rs)%sc(i-at(rs)%nbb)
        end if
        if (mass(ii).gt.5.0) then
          svr = svr + ln_ipol(ii)*atbvol(ii)
          if (do2.EQV..true.) svr2 = svr2 + ln_ipol2(ii)*atbvol(ii)
        end if 
      end do
!     ... by first normalizing with respect to the max. residue-SAV ...
      svr = svr/resvol(rs,2)
      svr2 = svr2/resvol(rs,2)
!     ... and by then multiplying with an appropriate reference model ...
!     ... "as is" for atoms
      if ((resname.eq.'CL-').OR.(resname.eq.'NA+').OR.&
 &        (resname.eq.'BR-').OR.(resname.eq.'K+ ').OR.&
 &        (resname.eq.'I- ').OR.(resname.eq.'CS+')) then
        tsv = tsv + svr*atsav(at(rs)%bb(1))
        if (do2.EQV..true.) then
          tsv2=tsv2 + svr2*(atsav(at(rs)%bb(1))+svte(at(rs)%bb(1)))
        end if
!     ... spherical assumption for small residues     
      else if ((resname.eq.'CH4').OR.(resname.eq.'SPC').OR.&
 &             (resname.eq.'FOR').OR.(resname.eq.'NH2').OR.&
 &             (resname.eq.'T3P').OR.(resname.eq.'T4P').OR.&
 &             (resname.eq.'T5P').OR.(resname.eq.'CH4').OR.&
 &             (resname.eq.'T4E').OR.&
 &             (resname.eq.'MOH').OR.(resname.eq.'NH4')) then
        sphrad = ((3./4.)*(1./PI)*resvol(rs,1))**(1./3.)
        sphvol = (4./3.)*PI*((sphrad+par_IMPSOLV(1))**3-sphrad**3)
        tsv = tsv + svr*sphvol
        if (do2.EQV..true.) then
          tsv2 = tsv2 + svr2*sphvol
        end if
!     ... and a cylindrical assumption for all larger residues
      else
        cylrad = resvol(rs,1)/(1.5*1.5*PI)
        cylvol = PI*(cylrad+par_IMPSOLV(1))*&
 &               ((1.5+par_IMPSOLV(1))**2.0)
        tsv = tsv + svr*cylvol
        if (do2.EQV..true.) then
          tsv2 = tsv2 + svr2*cylvol
        end if
      end if
    end do
  end do
!
end
!
!------------------------------------------------------------------------------
!
! a routine to set up the solvation groups, i.e., the atoms and weight factors
! contributing to each solvation unit along the chain
!
subroutine solvation_groups()
!
  use energies
  use sequen
  use polypep
  use aminos
  use iounit
  use molecule
  use fos
  use params
  use atoms
  use system
  use mpistuff
!
  implicit none
!
  integer i,j,k,t1,t2,rs,imol,dppdb,g,gg,freeunit,ii,jj,nps,fs
  character(3) resname,resname2
#ifdef ENABLE_MPI
  character(re_aux(10)) nod
#endif 
  character(200) fn
  character sty,sty2
  logical exists
  integer, ALLOCATABLE:: tmpmat(:,:),tmpmat3(:,:)
  RTYPE, ALLOCATABLE:: tmpmat2(:,:)
!  RTYPE fos_Tdep
!
  do imol=1,nmol
!
    do rs=rsmol(imol,1),rsmol(imol,2)
!
      resname = amino(seqtyp(rs))
!
      if (resname.eq.'UNK') at(rs)%nfosgrps = 0 ! may be overridden
!
      sty = seqpolty(rs)
      if (rs.gt.rsmol(imol,1)) then
        resname2 = amino(seqtyp(rs-1))
        sty2 = seqpolty(rs-1)
      end if
!
!     we'll start with the backbone-associated solvation groups only
!
!     for residues that represent a simple molecule (like ions, water) the treatment
!     is mostly straightforward

      if ((rsmol(imol,2)-rsmol(imol,1)).eq.0) then
!
        call helper_svgrps_single(imol,rs,resname,sty,sty2)
!
!
!     now let's turn our attention to C-terminal residues
!
      else if (rs.eq.rsmol(imol,2)) then
!
        call helper_svgrps_Cter(imol,rs,resname,resname2,sty,sty2)
!
!
!     now let's move on to N-terminal residues
! 
      else if (rs.eq.rsmol(imol,1)) then
!
        call helper_svgrps_Nter(imol,rs,resname,sty,sty2)
!
!
!     and finally take care of residues in the middle of the polymer chain
!
      else
!
        call helper_svgrps_middle(imol,rs,resname,resname2,sty,sty2)
!
      end if
!
!
!     the latter half of this function deals with the sidechain-associated
!     fos-groups
!
      call helper_svgrps_scpep(imol,rs,resname,sty,sty2)
      call helper_svgrps_scnuc(imol,rs,resname,sty,sty2)
!
!
    end do
!    
  end do
!
  allocate(tmpmat(n,2))
  tmpmat(:,:) = 0
  allocate(tmpmat2(n,4))
  tmpmat2(:,:) = 0.0
  allocate(tmpmat3(n,2))
  tmpmat3(:,:) = 0
!
  call read_fospatchfile(tmpmat,tmpmat2,nps)
!
! we get back a list of proposed FOS groups that satisfy the total weight constraint 
! now we need to check whether the boundaries overlap exactly with exisiting FOS group boundaries
! (tedious and inefficient but simple code)
 44 format(' Group ',i3,' in residue ',i6,' with DG_298 = ',g9.3,', DH = ',g9.3,', DCP = ',g9.3)
 45 format('    Atom  ',i7,': ',g10.4)
  if (nps.gt.0) then
    write(ilog,*)
    write(ilog,*) 'Input for FOS patches processed - Report follows'
    write(ilog,*)
    ii = 0
    do j=1,nps
      exists = .false.
      if (which_fosg(tmpmat(j,1)).le.0) then
        t1 = atmres(tmpmat(j,1))
        t2 = 0
        do i=1,ii
          if ((tmpmat3(i,1).eq.t1).AND.(tmpmat3(i,2).eq.t2)) then
            exists = .true.
            exit
          end if
        end do
        if (exists.EQV..false.) then
          ii = ii + 1
          tmpmat3(ii,1) = t1 ! atmres(tmpmat(j,1))
          tmpmat3(ii,2) = t2 ! which_fosg(tmpmat(j,1))
        end if
        cycle
      end if
      t1 = atmres(tmpmat(j,1))
      t2 = which_fosg(tmpmat(j,1))
      if (t2.gt.100) then
        t1 = t1 + 1
        t2 = t2 - 100
      end if
      do i=1,ii
        if ((tmpmat3(i,1).eq.t1).AND.(tmpmat3(i,2).eq.t2)) then
          exists = .true.
          exit
        end if
      end do
      if (exists.EQV..false.) then
        ii = ii + 1
        tmpmat3(ii,1) = t1 ! atmres(tmpmat(j,1))
        tmpmat3(ii,2) = t2 ! which_fosg(tmpmat(j,1))
      end if
    end do
    do i=1,ii
!     should be redundant 
      if (tmpmat3(i,2).gt.100) then
        tmpmat3(i,1) = tmpmat3(i,1)+1
        tmpmat3(i,2) = tmpmat3(i,2)-100
      end if
    end do
!   second round to catch residue-shift leftovers
    do j=1,nps
      exists = .false.
      if (which_fosg(tmpmat(j,1)).le.0) cycle
      t1 = atmres(tmpmat(j,1))
      t2 = which_fosg(tmpmat(j,1))
      if (t2.gt.100) t2 = 0
      do i=1,ii
        if ((tmpmat3(i,1).eq.t1).AND.(tmpmat3(i,2).eq.t2)) then
          exists = .true.
          exit
        end if
      end do
      if (exists.EQV..false.) then
        ii = ii + 1
        tmpmat3(ii,1) = t1 ! atmres(tmpmat(j,1))
        tmpmat3(ii,2) = t2 ! which_fosg(tmpmat(j,1))
      end if
    end do
    write(ilog,*) '--- Affected FOS Group Assignments before Patch ---'
    j = 0
    do i=1,ii
      if (tmpmat3(i,2).eq.0) cycle
      j = j + 1
      write(ilog,44) tmpmat3(i,2),tmpmat3(i,1),at(tmpmat3(i,1))%fosgrp(tmpmat3(i,2))%val(1:3)
      do k=1,at(tmpmat3(i,1))%fosgrp(tmpmat3(i,2))%nats
        write(ilog,45)  at(tmpmat3(i,1))%fosgrp(tmpmat3(i,2))%ats(k),at(tmpmat3(i,1))%fosgrp(tmpmat3(i,2))%wts(k)
      end do
    end do
    if (j.eq.0) then
      write(ilog,*) 'None'
      write(ilog,*) '---'
    else
      write(ilog,*) '---'
    end if
!
    call strlims(fospatchfile,t1,t2)
    do i=1,ii
      if (tmpmat3(i,2).eq.0) cycle
      do k=1,at(tmpmat3(i,1))%fosgrp(tmpmat3(i,2))%nats
        exists = .false.
        do j=1,nps
          if (at(tmpmat3(i,1))%fosgrp(tmpmat3(i,2))%ats(k).eq.tmpmat(j,1)) then
            exists = .true.
            exit
          end if
        end do
        if (exists.EQV..false.) then
          write(ilog,*) 'Fatal. Not all atoms in original FOS group (',tmpmat3(i,2),') in residue ',tmpmat3(i,1),&
 &' are accounted for by proposed FOS patch (file ',fospatchfile(t1:t2),'). Please correct input file.'
          call fexit()
        end if
      end do
    end do
!
!   now we can transfer the overrides
    fos_patched = .true.
    call helper_svgrps_patch(nps,ii,tmpmat,tmpmat2,tmpmat3)
!
    write(ilog,*) '--- Affected FOS Group Assignments after Patch  ---'
    i = 1
    do while (i.le.nps)
      t1 = tmpmat(i,2)
      fs = 0
      do while (tmpmat(i,2).eq.t1)
        if ((which_fosg(tmpmat(i,1)).gt.0).AND.(fs.le.0)) then
          fs = which_fosg(tmpmat(i,1))
          rs = atmres(tmpmat(i,1))
        end if
        i = i + 1
        if (i.gt.nps) exit
      end do
      if (fs.eq.0) then
        write(ilog,*) 'Fatal. This is a bug in solvation_groups(...). Please report this.'
        call fexit()
      end if
      if (fs.gt.100) then
        fs = fs - 100
        rs = rs + 1
      end if
      write(ilog,44) fs,rs, at(rs)%fosgrp(fs)%val(1:3)
      do k=1,at(rs)%fosgrp(fs)%nats
        write(ilog,45) at(rs)%fosgrp(fs)%ats(k),at(rs)%fosgrp(fs)%wts(k)
      end do
    end do
    do i=1,ii
      if (tmpmat3(i,1).gt.0) cycle
      if (tmpmat3(i,2).eq.0) cycle
      write(ilog,44) tmpmat3(i,2),-tmpmat3(i,1),at(-tmpmat3(i,1))%fosgrp(tmpmat3(i,2))%val(1:3)
      do k=1,at(-tmpmat3(i,1))%fosgrp(tmpmat3(i,2))%nats
        write(ilog,45)  at(-tmpmat3(i,1))%fosgrp(tmpmat3(i,2))%ats(k),at(-tmpmat3(i,1))%fosgrp(tmpmat3(i,2))%wts(k)
      end do
    end do
   write(ilog,*) '---'
!
  end if
!
  deallocate(tmpmat)
  deallocate(tmpmat2)
  deallocate(tmpmat3)
!  do rs=1,nseq
!    resname = amino(seqtyp(rs))
!    write(*,*) '###### ',rs,resname,':'
! 26    format(i3,1x,3(g10.4,1x),1x,i5)
! 27    format(3x,100(i3,2x))
! 28    format(3x,100(a3,2x))
! 29    format(3x,100(f5.2))
!    do j=1,at(rs)%nfosgrps
!      write(*,26) j,at(rs)%fosgrp(j)%val(1:3),at(rs)%fosgrp(j)%nats
!      write(*,27) (at(rs)%fosgrp(j)%ats(k),k=1,at(rs)%fosgrp(j)%nats)
!      write(*,28) (bio_code(b_type(at(rs)%fosgrp(j)%ats(k))),k=1,at(rs)%fosgrp(j)%nats)
!      write(*,27) (which_fosg(at(rs)%fosgrp(j)%ats(k)),k=1,at(rs)%fosgrp(j)%nats)
!      write(*,29) (at(rs)%fosgrp(j)%wts(k),k=1,at(rs)%fosgrp(j)%nats)
!    end do
!  end do
!  allocate(tmpmat2(2*sum(at(1:nseq)%nfosgrps),2))
! 555 format(100000(g11.4,1x))
! 556 format(100000(4x,a3,5x))
!  write(*,556) 'KEL',((amino(seqtyp(rs)),j=1,at(rs)%nfosgrps),rs=1,nseq)
!  do k=1,181
!    kelvin = 272.+ 1.0*k 
!    ii = 0
!    do rs=1,nseq
!      do j=1,at(rs)%nfosgrps
!         ii = ii + 1
!         tmpmat2(ii,1) = fos_Tdep(at(rs)%fosgrp(j)%val(1:3))
!      end do
!      if (at(rs)%nfosgrps.gt.2) then
!        ii = ii + 1
!        tmpmat2(ii,1) = sum(tmpmat2((ii-at(rs)%nfosgrps+1):(ii-1),1))
!      end if
!    end do
!    write(*,555) kelvin,tmpmat2(1:ii,1) !  /(kelvin*1.98717623d-3)
!  end do
!  call fexit()
!
  if (fos_report.EQV..true.) then
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then
      if (myrank.ne.0) return
      call int2str(myrank,nod,re_aux(10))
      fn = 'N_'//nod(1:re_aux(10))//'_FOS_GROUPS.vmd'
    else if (use_MPIAVG.EQV..true.) then
      fn = 'FOS_GROUPS.vmd'
      if (myrank.ne.0) return
    end if
#else
    fn = 'FOS_GROUPS.vmd'
#endif
    call strlims(fn,ii,jj)
    inquire(file=fn(ii:jj),exist=exists)
    if(exists) then
      dppdb = freeunit()
      open(unit=dppdb,file=fn(ii:jj),status='old')
      close(unit=dppdb,status='delete')
    end if
    dppdb=freeunit()
    open(unit=dppdb,file=fn(ii:jj),status='new')
!
    write(dppdb,*) 'color Display Background black'
    write(dppdb,*) 'material change opacity Transparent 0.52'
    write(dppdb,*) 'mol load pdb ',basename(1:bleng),'_START.pdb'
    write(dppdb,*) 'mol modcolor Type'
    write(dppdb,*) 'mol modstyle 0 0 Licorice'
!
    gg = 1
    do rs=1,nseq
!      write(*,*) amino(seqtyp(rs)),' with ',at(rs)%nfosgrps
      if (at(rs)%nfosgrps.gt.gg) gg = at(rs)%nfosgrps
!      do g=1,at(rs)%nfosgrps
!        damn = 0.0
!        do j=1,at(rs)%fosgrp(g)%nats
!          damn = damn + at(rs)%fosgrp(g)%wts(j)
!          write(*,*) bio_code(b_type(at(rs)%fosgrp(g)%ats(j)))
!        end do
!        write(*,*) damn
!      end do
!      damn = 0.0
!      do g=1,at(rs)%nfosgrps
!        damn = damn + at(rs)%fosgrp(g)%val(1)
!         write(*,*) g,at(rs)%fosgrp(g)%val(1)
!      end do
!      write(*,*) 'For ',damn
    end do
    do g=1,gg 
      write(dppdb,*) ' mol selection "index \'
      do rs=1,nseq
        if (at(rs)%nfosgrps.ge.g) then 
          do j=1,at(rs)%fosgrp(g)%nats
            write(dppdb,*) at(rs)%fosgrp(g)%ats(j)-1,'\'
          end do
        end if
      end do
      write(dppdb,*) '"'
      write(dppdb,*) ' mol color colorID ',g
      write(dppdb,*) ' mol rep VdW'
      write(dppdb,*) ' mol material Transparent'
      write(dppdb,*) ' mol addrep 0'
    end do
    close(unit=dppdb)
  end if
!
end
!
!-----------------------------------------------------------------------
!
! the following are asets of helper routines doing the actual setup work
! grouped by type of residue
!
!-----------------------------------------------------------------------
!
! this one deals with single-residue molecules (incl. free amino acids), excluding
! sidechains in general
!
subroutine helper_svgrps_single(imol,rs,resname,sty,sty2)
!
  use iounit
  use molecule
  use polypep
  use atoms
  use sequen
  use fos
  use zmatrix
  use system
!
  implicit none
!
  integer rs,imol,j,jj,hh,tt,shf,shf3,shf5,off,shfx1,shfl,shfl2
  character(3) resname
  character sty,sty2
!
  shf = 0
  shf3 = 0
  shf5 = 0
  shfx1 = 0
  if (ua_model.gt.0) then
    shf = 1
    shf3 = 3
    shf5 = 5
    if (ua_model.eq.2) then
      shfx1 = 1
    end if
  end if
!
! for residues that represent a simple molecule (like ions, water) the treatment
! is mostly straightforward

  if (resname.eq.'UNK') then
    write(ilog,*) "Warning. Solvation groups for free &
 &(nonpolymeric) and unsupported residue ",rs," are not set up (patch required)."
    return
  else if ((resname.eq.'NA+').OR.(resname.eq.'CL-').OR.&
 &    (resname.eq.'K+ ').OR.(resname.eq.'BR-').OR.&
 &    (resname.eq.'CS+').OR.(resname.eq.'I- ')) then
    at(rs)%nfosgrps = 1
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps))
    at(rs)%fosgrp(1)%val(1:3) = freesolv(rs,1:3)
    hh = 1
    allocate(at(rs)%fosgrp(1)%ats(hh))
    allocate(at(rs)%fosgrp(1)%wts(hh))
    at(rs)%fosgrp(1)%nats = hh
    at(rs)%fosgrp(1)%ats(1) = at(rs)%bb(1)
    at(rs)%fosgrp(1)%wts(1) = 1.0/(1.0*hh)
    which_fosg(at(rs)%bb(1)) = 1
!
  else if ((resname.eq.'SPC').OR.(resname.eq.'T3P'))  then
    at(rs)%nfosgrps = 1
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps))
    at(rs)%fosgrp(1)%val(1:3) = freesolv(rs,1:3)
    hh = at(rs)%nbb
    allocate(at(rs)%fosgrp(1)%ats(hh))
    allocate(at(rs)%fosgrp(1)%wts(hh))
    at(rs)%fosgrp(1)%nats = hh
    do j=1,at(rs)%nbb
      at(rs)%fosgrp(1)%ats(j) = at(rs)%bb(j)
      at(rs)%fosgrp(1)%wts(j) = 1.0/(1.0*hh)
      which_fosg(at(rs)%bb(j)) = 1
    end do
! for methanol treatment is analog. to SER/CYS (since the molecule is so small it doesn't
! really matter whether the CH3-unit is included)
  else if ((resname.eq.'MOH').OR.(resname.eq.'MSH')) then
    at(rs)%nfosgrps = 1
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps))
    at(rs)%fosgrp(1)%val(1:3) = freesolv(rs,1:3)
    hh = 2
    allocate(at(rs)%fosgrp(1)%ats(hh))
    allocate(at(rs)%fosgrp(1)%wts(hh))
    at(rs)%fosgrp(1)%nats = hh
    do j=1,min(2,at(rs)%nsc) ! MSH in ua_model == 2 has no SH
      at(rs)%fosgrp(1)%ats(j) = at(rs)%sc(j)
      at(rs)%fosgrp(1)%wts(j) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = 1
    end do
! note that for the acetate ion, the CH3 group is in sidechain array and hence not considered
! similarly, in 1MN and 2MN methyl hydrogens are sidechain as well
  else if ((resname.eq.'URE').OR.(resname.eq.'CH4').OR.&
 &         (resname.eq.'NH4').OR.(resname.eq.'AC-').OR.&
 &         (resname.eq.'GDN').OR.(resname.eq.'LCP').OR.&
 &         (resname.eq.'NO3').OR.(resname.eq.'1MN').OR.&
 &         (resname.eq.'2MN')) then
    at(rs)%nfosgrps = 1
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps))
    at(rs)%fosgrp(1)%val(1:3) = freesolv(rs,1:3)
    hh = at(rs)%nbb
    allocate(at(rs)%fosgrp(1)%ats(hh))
    allocate(at(rs)%fosgrp(1)%wts(hh))
    at(rs)%fosgrp(1)%nats = hh
    do j=1,at(rs)%nbb
      at(rs)%fosgrp(1)%ats(j) = at(rs)%bb(j)
      at(rs)%fosgrp(1)%wts(j) = 1.0/(1.0*hh)
      which_fosg(at(rs)%bb(j)) = 1
    end do
! isotropic ones
  else if ((resname.eq.'PRP').OR.(resname.eq.'TOL').OR.&
 &         (resname.eq.'IBU').OR.(resname.eq.'NBU').OR.&
 &         (resname.eq.'BEN').OR.(resname.eq.'NAP').OR.&
 &         (resname.eq.'O2 ')) then
    at(rs)%nfosgrps = 1
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps))
    at(rs)%fosgrp(1)%val(1:3) = freesolv(rs,1:3)
    hh = at(rs)%nbb + at(rs)%nsc
    allocate(at(rs)%fosgrp(1)%ats(hh))
    allocate(at(rs)%fosgrp(1)%wts(hh))
    at(rs)%fosgrp(1)%nats = hh
    do j=1,at(rs)%nbb+at(rs)%nsc
      if (j.le.at(rs)%nbb) then
        jj = at(rs)%bb(j)
      else
        jj = at(rs)%sc(j-at(rs)%nbb)
      end if
      at(rs)%fosgrp(1)%ats(j) = jj
      at(rs)%fosgrp(1)%wts(j) = 1.0/(1.0*hh)
      which_fosg(jj) = 1
    end do
! for the amides treatment is just like treatment for the backbone/sidechains
  else if (resname.eq.'FOA') then
    at(rs)%nfosgrps = 1
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps))
    at(rs)%fosgrp(1)%val(1:3) = freesolv(rs,1:3)
    hh = 5
    allocate(at(rs)%fosgrp(1)%ats(hh))
    allocate(at(rs)%fosgrp(1)%wts(hh))
    at(rs)%fosgrp(1)%nats = hh
    tt = 0
    do j=1,6-shf
      if ((ua_model.eq.0).AND.(j.eq.4)) cycle
      tt = tt + 1
      at(rs)%fosgrp(1)%ats(tt) = at(rs)%bb(j)
      at(rs)%fosgrp(1)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%bb(j)) = 1
    end do
  else if (resname.eq.'NMF') then
    at(rs)%nfosgrps = 1
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps))
    at(rs)%fosgrp(1)%val(1:3) = freesolv(rs,1:3)
    hh = 4
    allocate(at(rs)%fosgrp(1)%ats(hh))
    allocate(at(rs)%fosgrp(1)%wts(hh))
    at(rs)%fosgrp(1)%nats = hh
    tt = 0
    do j=1,5-shf
      if ((ua_model.eq.0).AND.(j.eq.4)) cycle
      tt = tt + 1
      at(rs)%fosgrp(1)%ats(tt) = at(rs)%bb(j)
      at(rs)%fosgrp(1)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%bb(j)) = 1
    end do
  else if (resname.eq.'NMA') then
    at(rs)%nfosgrps = 1
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps))
    at(rs)%fosgrp(1)%val(1:3) = freesolv(rs,1:3)
    hh = 4
    allocate(at(rs)%fosgrp(1)%ats(hh))
    allocate(at(rs)%fosgrp(1)%wts(hh))
    at(rs)%fosgrp(1)%nats = hh
    do j=1,4
      at(rs)%fosgrp(1)%ats(j) = at(rs)%bb(j)
      at(rs)%fosgrp(1)%wts(j) = 1.0/(1.0*hh)
      which_fosg(at(rs)%bb(j)) = 1
    end do
  else if ((resname.eq.'ACA').OR.(resname.eq.'PPA')) then
    at(rs)%nfosgrps = 1
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps))
    at(rs)%fosgrp(1)%val(1:3) = freesolv(rs,1:3)
    hh = 5
    allocate(at(rs)%fosgrp(1)%ats(hh))
    allocate(at(rs)%fosgrp(1)%wts(hh))
    at(rs)%fosgrp(1)%nats = hh
    do j=1,5
      at(rs)%fosgrp(1)%ats(j) = at(rs)%bb(j)
      at(rs)%fosgrp(1)%wts(j) = 1.0/(1.0*hh)
      which_fosg(at(rs)%bb(j)) = 1
    end do
! we'll follow the proline bb paradigm
  else if (resname.eq.'DMA') then
    at(rs)%nfosgrps = 1
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps))
    at(rs)%fosgrp(1)%val(1:3) = freesolv(rs,1:3)
    hh = 3
    allocate(at(rs)%fosgrp(1)%ats(hh))
    allocate(at(rs)%fosgrp(1)%wts(hh))
    at(rs)%fosgrp(1)%nats = hh
    do j=1,3
      at(rs)%fosgrp(1)%ats(j) = at(rs)%bb(j)
      at(rs)%fosgrp(1)%wts(j) = 1.0/(1.0*hh)
      which_fosg(at(rs)%bb(j)) = 1
    end do
!  for p-Cresol treatment is analog. to TYR
!  note that all decompoisiton schemes make sense for all 3 parameters due to additivity in the T-dependent extrapolation
  else if (resname.eq.'PCR') then
    at(rs)%nfosgrps = 2
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps))
!   the OH-part
    at(rs)%fosgrp(1)%val(1:3) = freesolv(rs,1:3) - FOS_PHE(1:3)
    hh = 2
    allocate(at(rs)%fosgrp(1)%ats(hh))
    allocate(at(rs)%fosgrp(1)%wts(hh))
    at(rs)%fosgrp(1)%nats = hh
    at(rs)%fosgrp(1)%ats(1) = at(rs)%bb(8)
    at(rs)%fosgrp(1)%wts(1) = 1.0/(1.0*hh)
    which_fosg(at(rs)%bb(8)) = 1
    at(rs)%fosgrp(1)%ats(2) = at(rs)%sc(4-shf3)
    at(rs)%fosgrp(1)%wts(2) = 1.0/(1.0*hh)
    which_fosg(at(rs)%sc(4-shf3)) = 1
!   the toluene-part (favorable by itsself!)
    at(rs)%fosgrp(2)%val(1:3) = FOS_PHE(1:3)
    hh = at(rs)%nbb + at(rs)%nsc - 2
    allocate(at(rs)%fosgrp(2)%ats(hh))
    allocate(at(rs)%fosgrp(2)%wts(hh))
    at(rs)%fosgrp(2)%nats = hh
    j = 0
    do tt=1,at(rs)%nbb
      if (tt.eq.8) cycle
      j = j + 1
      at(rs)%fosgrp(2)%ats(j) = at(rs)%bb(tt)
      at(rs)%fosgrp(2)%wts(j) = 1.0/(1.0*hh)
      which_fosg(at(rs)%bb(tt)) = 2
    end do
    if (ua_model.eq.0) then
      do tt=1,3
        j = j + 1 
        at(rs)%fosgrp(2)%ats(j) = at(rs)%sc(tt)
        at(rs)%fosgrp(2)%wts(j) = 1.0/(1.0*hh)
        which_fosg(at(rs)%sc(tt)) = 2
      end do
    end if
!  for EMT, just like MET
  else if (resname.eq.'EMT') then
    at(rs)%nfosgrps = 2
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps))
!   S
    at(rs)%fosgrp(1)%val(1:3) = freesolv(rs,1:3) - FOS_NBU(1:3) ! or FOS_NLE or FOS_ILE
    hh = 1
    allocate(at(rs)%fosgrp(1)%ats(hh))
    allocate(at(rs)%fosgrp(1)%wts(hh))
    at(rs)%fosgrp(1)%nats = hh
    at(rs)%fosgrp(1)%ats(1) = at(rs)%bb(2)
    at(rs)%fosgrp(1)%wts(1) = 1.0/(1.0*hh)
    which_fosg(at(rs)%bb(2)) = 1
!   the rest
    at(rs)%fosgrp(2)%val(1:3) = FOS_NBU(1:3)
    hh = at(rs)%nsc + at(rs)%nbb - 1
    allocate(at(rs)%fosgrp(2)%ats(hh))
    allocate(at(rs)%fosgrp(2)%wts(hh))
    at(rs)%fosgrp(2)%nats = hh
    tt = 0
    do j=1,at(rs)%nbb+at(rs)%nsc
      if (j.le.at(rs)%nbb) then
        jj = at(rs)%bb(j)
      else
        jj = at(rs)%sc(j-at(rs)%nbb)
      end if
      if (jj.eq.at(rs)%bb(2)) cycle
      tt = tt + 1
      at(rs)%fosgrp(2)%ats(tt) = jj
      at(rs)%fosgrp(2)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(jj) = 2
    end do
!  for EOH, just like THR
  else if (resname.eq.'EOH') then
    at(rs)%nfosgrps = 2
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps))
!   OH
    at(rs)%fosgrp(1)%val(1:3) = FOS_MOH(1:3) ! or FOS_SER
    hh = 2
    allocate(at(rs)%fosgrp(1)%ats(hh))
    allocate(at(rs)%fosgrp(1)%wts(hh))
    at(rs)%fosgrp(1)%nats = hh
    at(rs)%fosgrp(1)%ats(1) = at(rs)%sc(4-shf3)
    at(rs)%fosgrp(1)%wts(1) = 1.0/(1.0*hh)
    which_fosg(at(rs)%sc(4-shf3)) = 1
    at(rs)%fosgrp(1)%ats(2) = at(rs)%bb(3)
    at(rs)%fosgrp(1)%wts(2) = 1.0/(1.0*hh)
    which_fosg(at(rs)%bb(3)) = 1
!   CH3
    at(rs)%fosgrp(2)%val(1:3) = freesolv(rs,1:3) - FOS_MOH(1:3)
    hh = 4 - shf3
    allocate(at(rs)%fosgrp(2)%ats(hh))
    allocate(at(rs)%fosgrp(2)%wts(hh))
    at(rs)%fosgrp(2)%nats = hh
    at(rs)%fosgrp(2)%nats = hh
    at(rs)%fosgrp(2)%ats(1) = at(rs)%bb(1)
    at(rs)%fosgrp(2)%wts(1) = 1.0/(1.0*hh)
    which_fosg(at(rs)%bb(1)) = 2
    if (ua_model.eq.0) then
      do j=1,3
        at(rs)%fosgrp(2)%ats(j+1) = at(rs)%sc(j)
        at(rs)%fosgrp(2)%wts(j+1) = 1.0/(1.0*hh)
        which_fosg(at(rs)%sc(j)) = 2
      end do
    end if
!  for MIN, just like TRP
  else if (resname.eq.'MIN') then
    at(rs)%nfosgrps = 2
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps))
!   NH
    at(rs)%fosgrp(1)%val(1:3) = freesolv(rs,1:3) - FOS_NAP(1:3)
    hh = 2
    allocate(at(rs)%fosgrp(1)%ats(hh))
    allocate(at(rs)%fosgrp(1)%wts(hh))
    at(rs)%fosgrp(1)%nats = hh
    at(rs)%fosgrp(1)%ats(1) = at(rs)%bb(5)
    at(rs)%fosgrp(1)%wts(1) = 1.0/(1.0*hh)
    which_fosg(at(rs)%bb(5)) = 1
    at(rs)%fosgrp(1)%ats(2) = at(rs)%bb(12-shfx1)
    at(rs)%fosgrp(1)%wts(2) = 1.0/(1.0*hh)
    which_fosg(at(rs)%bb(12-shfx1)) = 1
!   rest
    at(rs)%fosgrp(2)%val(1:3) = FOS_NAP(1:3)
    hh = at(rs)%nsc + at(rs)%nbb - 2
    allocate(at(rs)%fosgrp(2)%ats(hh))
    allocate(at(rs)%fosgrp(2)%wts(hh))
    at(rs)%fosgrp(2)%nats = hh
    tt = 0
    do j=1,at(rs)%nbb+at(rs)%nsc
      if (j.le.at(rs)%nbb) then
        jj = at(rs)%bb(j)
      else
        jj = at(rs)%sc(j-at(rs)%nbb)
      end if
      if ((jj.eq.at(rs)%bb(5)).OR.(jj.eq.at(rs)%bb(12-shfx1))) cycle
      tt = tt + 1
      at(rs)%fosgrp(2)%ats(tt) = jj
      at(rs)%fosgrp(2)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(jj) = 2
    end do
!  for IMX, just like HIX
  else if ((resname.eq.'IME').OR.(resname.eq.'IMD')) then
    at(rs)%nfosgrps = 1
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps))
!   HNCN
    hh = 4
    allocate(at(rs)%fosgrp(1)%ats(hh))
    allocate(at(rs)%fosgrp(1)%wts(hh))
    at(rs)%fosgrp(1)%val(1:3) = freesolv(rs,1:3)
    at(rs)%fosgrp(1)%nats = hh
    at(rs)%fosgrp(1)%ats(1) = at(rs)%bb(3)
    at(rs)%fosgrp(1)%wts(1) = 1.0/(1.0*hh)
    which_fosg(at(rs)%bb(3)) = 1
    at(rs)%fosgrp(1)%ats(2) = at(rs)%bb(5)
    at(rs)%fosgrp(1)%wts(2) = 1.0/(1.0*hh)
    which_fosg(at(rs)%bb(5)) = 1
    at(rs)%fosgrp(1)%ats(3) = at(rs)%bb(6)
    at(rs)%fosgrp(1)%wts(3) = 1.0/(1.0*hh)
    which_fosg(at(rs)%bb(6)) = 1
    off = 0
    if ((resname.eq.'IME').AND.(ua_model.lt.2)) off = 2
    at(rs)%fosgrp(1)%ats(4) = at(rs)%bb(7+off)
    at(rs)%fosgrp(1)%wts(4) = 1.0/(1.0*hh)
    which_fosg(at(rs)%bb(7+off)) = 1
!   for ADE/PUR, similar to DPA
  else if ((resname.eq.'ADE').OR.(resname.eq.'PUR')) then
    shfl = 1
    if (resname.eq.'ADE') shfl = 0
    at(rs)%nfosgrps = 3 - shfl
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps))
!   the polar NH
    hh = 2
    allocate(at(rs)%fosgrp(1)%ats(hh))
    allocate(at(rs)%fosgrp(1)%wts(hh))
    at(rs)%fosgrp(1)%nats = hh
    at(rs)%fosgrp(1)%val(1:3) = FOS_ADE(1:3) - FOS_XPA(1:3)
    at(rs)%fosgrp(1)%ats(1) = at(rs)%bb(1)
    at(rs)%fosgrp(1)%wts(1) = 1.0/(1.0*hh)
    which_fosg(at(rs)%bb(1)) = 1
    at(rs)%fosgrp(1)%ats(2) = at(rs)%bb(11-shfl)
    at(rs)%fosgrp(1)%wts(2) = 1.0/(1.0*hh)
    which_fosg(at(rs)%bb(11-shfl)) = 1
!   the main ring
    shfl2 = 2 + shfl
    if (ua_model.eq.2) shfl2 = 0
    hh = 8 + shfl2
    allocate(at(rs)%fosgrp(2)%ats(hh))
    allocate(at(rs)%fosgrp(2)%wts(hh))
    at(rs)%fosgrp(2)%nats = hh
    at(rs)%fosgrp(2)%val(1:3) = FOS_PUR(1:3)
    tt = 0
    do j=2,at(rs)%nbb+2*(shfl-1)
      if ((j.eq.(11-shfl)).OR.(j.eq.10)) cycle
      tt = tt + 1
      at(rs)%fosgrp(2)%ats(tt) = at(rs)%bb(j)
      at(rs)%fosgrp(2)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%bb(j)) = 2
    end do
    if (resname.eq.'ADE') then
!     the NH2-group
      hh = 3
      allocate(at(rs)%fosgrp(3)%ats(hh))
      allocate(at(rs)%fosgrp(3)%wts(hh))
      at(rs)%fosgrp(3)%nats = hh
      at(rs)%fosgrp(3)%val(1:3) = FOS_XPA(1:3) - FOS_PUR(1:3)
      tt = 0
      do j=at(rs)%bb(10),at(rs)%bb(at(rs)%nbb)
        if ((iz(1,j).eq.at(rs)%bb(10)).OR.(j.eq.at(rs)%bb(10))) then
          tt = tt + 1
          at(rs)%fosgrp(3)%ats(tt) = j
          at(rs)%fosgrp(3)%wts(tt) = 1.0/(1.0*hh)
          which_fosg(j) = 3
        end if
      end do
    end if
!   for URA, analogous to DPU
  else if (resname.eq.'URA') then
    shfl2 = 0
    if (ua_model.eq.2) shfl2 = 1
    at(rs)%nfosgrps = 2
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps))
!   the polar ring component
    hh = 8
    allocate(at(rs)%fosgrp(1)%ats(hh))
    allocate(at(rs)%fosgrp(1)%wts(hh))
    at(rs)%fosgrp(1)%nats = hh
    at(rs)%fosgrp(1)%val(1:3) = FOS_URA(1:3) - FOS_BEN(1:3)
    tt = 0
    do j=1,at(rs)%nbb
      if ((j.eq.3).OR.(j.eq.6)) cycle
      if ((j.eq.(11-shfl2)).OR.(j.le.9)) then
        tt = tt + 1
        at(rs)%fosgrp(1)%ats(tt) = at(rs)%bb(j)
        at(rs)%fosgrp(1)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%bb(j)) = 1
      end if
    end do
!   the apolar component
    hh = 4 - 2*shfl2
    allocate(at(rs)%fosgrp(2)%ats(hh))
    allocate(at(rs)%fosgrp(2)%wts(hh))
    at(rs)%fosgrp(2)%nats = hh
    at(rs)%fosgrp(2)%val(1:3) = FOS_BEN(1:3)
    tt = 0
    do j=3,at(rs)%nbb
      if ((j.eq.3).OR.(j.eq.6).OR.(j.eq.12).OR.((j.eq.10).AND.(shfl2.eq.0))) then
        tt = tt + 1
        at(rs)%fosgrp(2)%ats(tt) = at(rs)%bb(j)
        at(rs)%fosgrp(2)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%bb(j)) = 2
      end if
    end do
!   for THY, analogous to DPT
  else if (resname.eq.'THY') then
    shfl2 = 0
    if (ua_model.eq.2) shfl2 = 1
    at(rs)%nfosgrps = 3
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps))
!   the polar ring component
    hh = 8
    allocate(at(rs)%fosgrp(1)%ats(hh))
    allocate(at(rs)%fosgrp(1)%wts(hh))
    at(rs)%fosgrp(1)%nats = hh
    at(rs)%fosgrp(1)%val(1:3) = FOS_URA(1:3) - FOS_BEN(1:3)
    tt = 0
    do j=1,at(rs)%nbb
      if ((j.eq.3).OR.(j.eq.6).OR.(j.eq.8)) cycle
      if ((j.eq.(12-shfl2)).OR.(j.le.10)) then
        tt = tt + 1
        at(rs)%fosgrp(1)%ats(tt) = at(rs)%bb(j)
        at(rs)%fosgrp(1)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%bb(j)) = 1
      end if
    end do
!   the apolar component
    hh = 3 - shfl2
    allocate(at(rs)%fosgrp(2)%ats(hh))
    allocate(at(rs)%fosgrp(2)%wts(hh))
    at(rs)%fosgrp(2)%nats = hh
    at(rs)%fosgrp(2)%val(1:3) = FOS_BEN(1:3)
    tt = 0
    do j=3,at(rs)%nbb
      if ((j.eq.3).OR.(j.eq.6).OR.((j.eq.11).AND.(shfl2.eq.0))) then
        tt = tt + 1
        at(rs)%fosgrp(2)%ats(tt) = at(rs)%bb(j)
        at(rs)%fosgrp(2)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%bb(j)) = 2
      end if
    end do
!   the methyl group
    hh = 4 - shf3
    allocate(at(rs)%fosgrp(3)%ats(hh))
    allocate(at(rs)%fosgrp(3)%wts(hh))
    at(rs)%fosgrp(3)%nats = hh
    at(rs)%fosgrp(3)%val(1:3) = FOS_THY(1:3) - FOS_URA(1:3)
    tt = 1
    at(rs)%fosgrp(3)%ats(tt) = at(rs)%bb(8)
    at(rs)%fosgrp(3)%wts(tt) = 1.0/(1.0*hh)
    which_fosg(at(rs)%bb(8)) = 3
    do j=1,at(rs)%nsc
      tt = tt + 1
      at(rs)%fosgrp(3)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(3)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = 3
    end do
!   for CYT, analogous to DPC
  else if (resname.eq.'CYT') then
    shfl2 = 0
    if (ua_model.eq.2) shfl2 = 1
    at(rs)%nfosgrps = 3
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps))
!   the polar ring component
    hh = 5
    allocate(at(rs)%fosgrp(1)%ats(hh))
    allocate(at(rs)%fosgrp(1)%wts(hh))
    at(rs)%fosgrp(1)%nats = hh
    at(rs)%fosgrp(1)%val(1:3) = FOS_URA(1:3) - FOS_BEN(1:3)
    tt = 0
    do j=1,at(rs)%nbb
      if ((j.le.2).OR.(j.eq.4).OR.(j.eq.5).OR.(j.eq.9)) then
        tt = tt + 1
        at(rs)%fosgrp(1)%ats(tt) = at(rs)%bb(j)
        at(rs)%fosgrp(1)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%bb(j)) = 1
      end if
    end do
!   the apolar component
    hh = 5 - 2*shfl2
    allocate(at(rs)%fosgrp(2)%ats(hh))
    allocate(at(rs)%fosgrp(2)%wts(hh))
    at(rs)%fosgrp(2)%nats = hh
    at(rs)%fosgrp(2)%val(1:3) = FOS_BEN(1:3)
    tt = 0
    do j=3,at(rs)%nbb
      if ((j.eq.3).OR.(j.eq.6).OR.(j.eq.7).OR.(((j.eq.10).OR.(j.eq.11)).AND.(shfl2.eq.0))) then
        tt = tt + 1
        at(rs)%fosgrp(2)%ats(tt) = at(rs)%bb(j)
        at(rs)%fosgrp(2)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%bb(j)) = 2
      end if
    end do
!   the aniline group
    hh = 3
    allocate(at(rs)%fosgrp(3)%ats(hh))
    allocate(at(rs)%fosgrp(3)%wts(hh))
    at(rs)%fosgrp(3)%nats = hh
    at(rs)%fosgrp(3)%val(1:3) = FOS_CYT(1:3) - FOS_URA(1:3)
    tt = 0
    do j=8,at(rs)%nbb
      if ((j.eq.8).OR.(j.ge.(at(rs)%nbb-1))) then
        tt = tt + 1
        at(rs)%fosgrp(3)%ats(tt) = at(rs)%bb(j)
        at(rs)%fosgrp(3)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%bb(j)) = 3
      end if
    end do
!   for GUA, analogous to DPG
  else if (resname.eq.'GUA') then
    shfl2 = 0
    if (ua_model.eq.2) shfl2 = 1
    at(rs)%nfosgrps = 4
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps))
!   the 9N-H
    hh = 2
    allocate(at(rs)%fosgrp(1)%ats(hh))
    allocate(at(rs)%fosgrp(1)%wts(hh))
    at(rs)%fosgrp(1)%nats = hh
    at(rs)%fosgrp(1)%val(1:3) = FOS_GUA(1:3) - FOS_XPG(1:3)
    tt = 0
    do j=1,at(rs)%nbb
      if ((j.eq.1).OR.(j.eq.12)) then
        tt = tt + 1
        at(rs)%fosgrp(1)%ats(tt) = at(rs)%bb(j)
        at(rs)%fosgrp(1)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%bb(j)) = 1
      end if
    end do
!   the adenine-like ring component
    hh = 7 - shfl2
    allocate(at(rs)%fosgrp(2)%ats(hh))
    allocate(at(rs)%fosgrp(2)%wts(hh))
    at(rs)%fosgrp(2)%nats = hh
    at(rs)%fosgrp(2)%val(1:3) = FOS_PUR(1:3)
    tt = 0
    do j=2,13-shfl2
      if ((j.eq.8).OR.(j.eq.9).OR.(j.eq.10).OR.(j.eq.11).OR.(j.eq.12)) cycle
      tt = tt + 1
      at(rs)%fosgrp(2)%ats(tt) = at(rs)%bb(j)
      at(rs)%fosgrp(2)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%bb(j)) = 2
    end do
!   the H-bonding ring component
    hh = 4
    allocate(at(rs)%fosgrp(3)%ats(hh))
    allocate(at(rs)%fosgrp(3)%wts(hh))
    at(rs)%fosgrp(3)%nats = hh
    at(rs)%fosgrp(3)%val(1:3) = FOS_XPG(1:3) - FOS_XPA(1:3)
    tt = 0
    do j=3,at(rs)%nbb
      if ((j.eq.8).OR.(j.eq.9).OR.(j.eq.11).OR.(j.eq.(14-shfl2))) then
        tt = tt + 1
        at(rs)%fosgrp(3)%ats(tt) = at(rs)%bb(j)
        at(rs)%fosgrp(3)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%bb(j)) = 3
      end if
    end do
!   the aniline group
    hh = 3
    allocate(at(rs)%fosgrp(4)%ats(hh))
    allocate(at(rs)%fosgrp(4)%wts(hh))
    at(rs)%fosgrp(4)%nats = hh
    at(rs)%fosgrp(4)%val(1:3) = FOS_XPA(1:3) - FOS_PUR(1:3)
    tt = 0
    do j=10,at(rs)%nbb
      if ((j.eq.10).OR.(j.ge.(at(rs)%nbb-1))) then
        tt = tt + 1
        at(rs)%fosgrp(4)%ats(tt) = at(rs)%bb(j)
        at(rs)%fosgrp(4)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%bb(j)) = 4
      end if
    end do
! backbone of free amino acids
  else if ((seqflag(rs).eq.2).OR.(seqflag(rs).eq.5).OR.(seqflag(rs).eq.8)) then
!
!   we'll overallocate and then just use however many we need
    allocate(at(rs)%fosgrp(5))
    at(rs)%nfosgrps = 0
!
    if (moltermid(molofrs(rs),2).eq.1) then
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PEP_CCT(1:3)
      hh = 3
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      tt = 0
      do j=3,5
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%bb(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%bb(j)) = at(rs)%nfosgrps
      end do
    else
      write(ilog,*) 'Fatal. Encountered unsupported C-terminus&
 & in setup_fosgrps() (offending mode is ',moltermid(molofrs(rs),2),&
 &').'
      call fexit()
    end if
!
    if (moltermid(molofrs(rs),1).eq.1) then
!
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      if (seqflag(rs).eq.5) then
        at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PEP_PROCNT(1:3)
        hh = 3
      else
        at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PEP_CNT(1:3)
        hh = 4
      end if
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      tt = 0
      if (seqflag(rs).eq.5) then
        do j=1,7
          if ((j.ge.2).AND.(j.le.5)) cycle
          tt = tt + 1
          at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt)= at(rs)%bb(j)
          at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt)= 1.0/(1.0*hh)
          which_fosg(at(rs)%bb(j)) = at(rs)%nfosgrps
        end do
      else
        do j=1,8
          if ((j.ge.2).AND.(j.le.5)) cycle
          tt = tt + 1
          at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt)= at(rs)%bb(j)
          at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt)= 1.0/(1.0*hh)
          which_fosg(at(rs)%bb(j)) = at(rs)%nfosgrps
        end do
      end if
!
    else if (moltermid(molofrs(rs),1).eq.2) then
!
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      if (seqflag(rs).eq.5) then
        at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PEP_PROUNT(1:3)
        hh = 2
      else
        at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PEP_UNT(1:3)
        hh = 3
      end if
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      tt = 0
      if (seqflag(rs).eq.5) then
        do j=1,6
          if ((j.ge.2).AND.(j.le.5)) cycle
          tt = tt + 1
          at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt)= at(rs)%bb(j)
          at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt)= 1.0/(1.0*hh)
          which_fosg(at(rs)%bb(j)) = at(rs)%nfosgrps
        end do
      else
        do j=1,7
          if ((j.ge.2).AND.(j.le.5)) cycle
          tt = tt + 1
          at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt)= at(rs)%bb(j)
          at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt)= 1.0/(1.0*hh)
          which_fosg(at(rs)%bb(j)) = at(rs)%nfosgrps
        end do
      end if
!
    else
      write(ilog,*) 'Fatal. Encountered unsupported N-term. in helpe&
 &r_svgrps_single(...) (offending mode is ',moltermid(imol,1),').'
      call fexit()
    end if
  else
    write(ilog,*) 'Fatal. Implicit solvent does not support re&
 &sidue ',resname,' (free) at the moment. Check back later.'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
! this one deals with C-terminal residues (except sidechains)
!
subroutine helper_svgrps_Cter(imol,rs,resname,resname2,sty,sty2)
!
  use iounit
  use molecule
  use polypep
  use atoms
  use sequen
  use fos
  use system
!
  implicit none
!
  integer rs,imol,j,hh,tt,shf,shf2,shf3,shf6,shf7
  character(3) resname,resname2
  character sty,sty2
!
  shf = 0
  shf2 = 0
  shf3 = 0
  shf6 = 0
  shf7 = 0
  if (ua_model.gt.0) then
    shf = 1
    shf2 = 2
    shf3 = 3
    shf6 = 6
    shf7 = 7
  end if
!
  if ((sty.eq.'P').OR.(seqtyp(rs).eq.29).OR.(seqtyp(rs).eq.30)) then
!
    if ((ci(rs-1).le.0).OR.(oi(rs-1).le.0)) then
      if (resname2.eq.'UNK') then
        write(ilog,*) 'Warning. Unusual polymer linkage due to unsupported residue leads to elimination of &
 &backbone solvation group for residue ',rs,' in helper_svgrps_Cter(...).'
        at(rs)%nfosgrps = 0
        return
      else
        write(ilog,*) 'Fatal. Reference atoms missing for assumed polypeptide linkage &
 &for the backbone solvation group for residue ',rs,' in helper_svgrps_middle(...). This is a bug.'
        call fexit()
      end if
    end if
!
    if ((sty2.ne.'P').AND.(seqflag(rs-1).ne.10)) then
      write(ilog,*) 'Fatal. Implicit solvent does not support this typ&
 &e of polymer (linkage) or these polypeptide caps at the moment.'
      call fexit()
    end if
!
!   for residues that are part of a polypeptide chain, pursue the following strategy
!   first the backbone, rs-1,rs peptide unit "belongs" to residue rs
!   the relevant atoms are given by the polar CONH skeleton, C-alpha solvation is deemed to be negligible
!   this approach more or less automatically takes care of capped termini (ACE,NME,FOR,NH2)
!
!   since the which_fosg-array normally assumes the fos-group is part of that residue we'll use
!   an offset to indicate atoms which belong to the group #x-100 in residue rs+1
!
!   we'll overallocate and then just use however many we need
    allocate(at(rs)%fosgrp(5))
    at(rs)%nfosgrps = 0
!
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    if (resname.eq.'NH2') then
      hh = 5
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
!
!     first add CO on the previous residue
!     note we use a shift-by-100 trick to indicate that ci(rs-1) and oi(rs-1)
!     are not part of residue rs-1's fos-group list
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1)= ci(rs-1)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1)= 1.0/(1.0*hh)
      which_fosg(ci(rs-1)) = at(rs)%nfosgrps + 100
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2)= oi(rs-1)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2)= 1.0/(1.0*hh)
      which_fosg(oi(rs-1)) = at(rs)%nfosgrps + 100
!     set the right ref. value
      if (resname2.eq.'FOR') then
        at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_FOA(1:3)
      else if (resname2.eq.'ACE') then
        at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_ACA(1:3)
      else
        at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PEP_BB_NH2(1:3)
      end if
!     add the NH2-atoms
      do j=1,3
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(j+2) = at(rs)%bb(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(j+2) = 1.0/(1.0*hh)
        which_fosg(at(rs)%bb(j)) = at(rs)%nfosgrps
      end do
!
    else if ((seqflag(rs).eq.2).OR.(seqflag(rs).eq.5).OR.(seqflag(rs).eq.8).OR.(resname.eq.'NME')) then
      if ((seqflag(rs).eq.5).OR.(hni(rs).le.0)) then
        hh = 3
      else
        hh = 4
      end if
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
!
!     first add CO on the previous residue
!     note we use a shift-by-100 trick to indicate that ci(rs-1) and oi(rs-1)
!     are not part of residue rs-1's fos-group list
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1)= ci(rs-1)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1)= 1.0/(1.0*hh)
      which_fosg(ci(rs-1)) = at(rs)%nfosgrps + 100
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2)= oi(rs-1)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2)= 1.0/(1.0*hh)
      which_fosg(oi(rs-1)) = at(rs)%nfosgrps + 100
      if (resname.eq.'NME') then
!       set the right ref. value
        if (resname2.eq.'FOR') then
          at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NMF(1:3)
        else if (resname2.eq.'ACE') then
          at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NMA(1:3)
        else
          at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PEP_BB(1:3)
        end if
!       add the NH-atoms
        tt = 2
        do j=1,3
          if (j.eq.2) cycle
          tt = tt + 1
          at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt)= at(rs)%bb(j)
          at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt)= 1.0/(1.0*hh)
          which_fosg(at(rs)%bb(j)) = at(rs)%nfosgrps
        end do
      else
!       set the right ref. value
        if ((seqflag(rs).eq.5).OR.(hni(rs).le.0)) then
          if (resname2.eq.'FOR') then
            at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PEP_PROBB_FOR(1:3)
          else
            at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PEP_PROBB(1:3)
          end if
        else
          if (resname2.eq.'FOR') then
            at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PEP_BB_FOR(1:3)
          else
            at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PEP_BB(1:3)
          end if
        end if
!       add the NH-atoms
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(3)= ni(rs)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(3)= 1.0/(1.0*hh)
        which_fosg(ni(rs)) = at(rs)%nfosgrps
        if ((seqflag(rs).ne.5).AND.(hni(rs).gt.0)) then
          at(rs)%fosgrp(at(rs)%nfosgrps)%ats(4)= hni(rs)
          at(rs)%fosgrp(at(rs)%nfosgrps)%wts(4)= 1.0/(1.0*hh)
          which_fosg(hni(rs)) = at(rs)%nfosgrps
        end if
      end if
    else if (seqtyp(rs).eq.26) then
      write(ilog,*) 'Warning. Even though unsupported residue is recognized as a C-terminal polypeptide cap, &
 &terminal solvation group(s) for residue ',rs,' are not set up (patch required).'
      return
    else
      write(ilog,*) 'Fatal. Implicit solvent does not support &
 &C-terminal residue ',resname,' at the moment. Check back later.'
      call fexit()
    end if
!
    if (sty.eq.'P') then
      if (resname.eq.'UNK') then
        write(ilog,*) 'Warning. Even though unsupported residue is recognized as polypeptide, terminal solvation &
 &group(s) for residue ',rs,' are not set up (patch required).'
        return
      end if
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      if (moltermid(imol,2).eq.1) then
        hh = 3
        allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
        allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
        at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
        at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PEP_CCT(1:3)
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1)= ci(rs)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1)= 1.0/(1.0*hh)
        which_fosg(ci(rs)) = at(rs)%nfosgrps
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2)= oi(rs)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2)= 1.0/(1.0*hh)
        which_fosg(oi(rs)) = at(rs)%nfosgrps
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(3)= at(rs)%bb(5)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(3)= 1.0/(1.0*hh)
        which_fosg(at(rs)%bb(5)) = at(rs)%nfosgrps
      else
        write(ilog,*) 'Fatal. C-terminus type ',&
 &moltermid(imol,2),'(mol. ',imol,') is not yet supported by implici&
 &t solvent model.'
        call fexit()
      end if
    end if
!
! for residues which are part of a polynucleotide chain we have the following problem:
! the XYphosphate (XYPO4-) and the ribose are esterified and hence both deprived of at least
! one H-bond donor at the ester site:
! similar to glycosylation (see below) this makes identification of clean model compounds
! difficult, as we practically have stunted model compounds only. the decomposition approach 
! followed here is to assume two different phosphate compounds (mono- and dimethylester) as well
! as interpolations between ribose and 2-methyl-THF for the sugar
! 
  else if (sty.eq.'N') then
!
    if (resname.eq.'UNK') then
      write(ilog,*) "Warning. Even though unsupported residue is recognized as nucleic acid, solvation &
 &groups for 3'-terminal residue ",rs," are not set up (patch required)."
      return
    end if
!
!   we'll overallocate and then just use however many we need
    allocate(at(rs)%fosgrp(7))
    at(rs)%nfosgrps = 0
!
!   first the phospho-diester (CCPO4-): we'll just use the POO skeleton in the middle
!   since the estered oxygens cannot carry charge (but overlap large anyway)
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 3
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_PO4D(1:3)
    tt = 0
    do j=2,5
      if (j.eq.3) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%bb(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%bb(j)) = at(rs)%nfosgrps
    end do
!
!   now the sugar: THF + decorating hydroxyls are really the underlying model
!   compounds for us: we'll start with the 3'-OH since it's common to all
    if (moltermid(imol,2).ne.1) then
      write(ilog,*) 'Fatal. C-terminus type ',&
 &moltermid(imol,2),'(mol. ',imol,') is not yet supported by implici&
 &it solvent model.'
      call fexit()
    end if
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 2
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_RIBO1(1:3) - FOS_NUC_MTHF(1:3)
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1) = at(rs)%bb(9)
    at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1) = 1.0/(1.0*hh)
    which_fosg(at(rs)%bb(9)) = at(rs)%nfosgrps
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2) = at(rs)%bb(at(rs)%nbb)
    at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2) = 1.0/(1.0*hh)
    which_fosg(at(rs)%bb(at(rs)%nbb)) = at(rs)%nfosgrps
!   now we add the rest: the ring core and 0-2 hydroxyl groups
    if ((resname.eq.'DPC').OR.(resname.eq.'DPU').OR.&
 &      (resname.eq.'DPT').OR.(resname.eq.'DPA').OR.&
 &      (resname.eq.'DPG')) then
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      hh = 13 - shf7 ! in UA, 2xH5*, H4*, H3*, 2xH2*, H1* are missing
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_MTHF(1:3)
      tt = 0
      do j=6,at(rs)%nbb-1
        if (j.eq.9) cycle
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%bb(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%bb(j)) = at(rs)%nfosgrps
      end do
      do j=1,6-shf3
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
      end do
    else if ((resname.eq.'RPC').OR.(resname.eq.'RPU').OR.&
 &           (resname.eq.'RPT').OR.(resname.eq.'RPA').OR.&
 &           (resname.eq.'RPG')) then
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      hh = 12 - shf6 ! in UA, 2xH5*, H4*, H3*, H2*, H1* are missing
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_MTHF(1:3)
      tt = 0
      do j=6,at(rs)%nbb-1
        if (j.eq.9) cycle
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%bb(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%bb(j)) = at(rs)%nfosgrps
      end do
      do j=1,6-shf3
        if (j.eq.4) cycle
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
      end do
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      hh = 2
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_RIBO2(1:3) - FOS_NUC_RIBO1(1:3)
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1) = at(rs)%sc(4)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(4)) = at(rs)%nfosgrps
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2) = at(rs)%sc(7-shf2)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(7-shf2)) = at(rs)%nfosgrps
    else if (resname.eq.'D5P') then
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      hh = 13 - shf7 ! in UA, 2xH5*, H4*, H3*, 2xH2*, H1* are missing
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_MTHF(1:3)
      tt = 0
      do j=6,at(rs)%nbb-1
        if (j.eq.9) cycle
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%bb(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%bb(j)) = at(rs)%nfosgrps
      end do
      do j=1,6-shf3
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
      end do
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      hh = 2
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_RIBO2(1:3) - FOS_NUC_RIBO1(1:3)

      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1) = at(rs)%sc(8-shf3)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(8-shf3)) = at(rs)%nfosgrps
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2) = at(rs)%sc(7-shf3)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(7-shf3)) = at(rs)%nfosgrps
    else if (resname.eq.'R5P') then
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      hh = 12 - shf6
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_MTHF(1:3)
      tt = 0
      do j=6,at(rs)%nbb-1
        if (j.eq.9) cycle
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%bb(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%bb(j)) = at(rs)%nfosgrps
      end do
      do j=1,6-shf3
        if (j.eq.4) cycle
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
      end do
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      hh = 2
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_RIBO2(1:3) - FOS_NUC_RIBO1(1:3)
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1) = at(rs)%sc(4)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(4)) = at(rs)%nfosgrps
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2) = at(rs)%sc(7-shf2)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(7-shf2)) = at(rs)%nfosgrps
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      hh = 2
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_RIBO3(1:3) - FOS_NUC_RIBO2(1:3)
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1) = at(rs)%sc(8-shf2)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(8-shf2)) = at(rs)%nfosgrps
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2) = at(rs)%sc(9-shf2)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(9-shf2)) = at(rs)%nfosgrps
    end if
!
  else if (resname.eq.'UNK') then
    write(ilog,*) 'Warning. C-terminal cap residue of unrecognized polymer type will &
 &have no solvation groups set up (patch required for residue ',rs,').'
    return
!
  else
!
    write(ilog,*) 'Fatal. Implicit solvent does not support this typ&
 &e of polymer or these polypeptide caps at the moment.'
    call fexit()
!
  end if
!
end
!
!-----------------------------------------------------------------------
!
! surprise, surprise, this one deals with N-terminal residues (except sidechains)
!
subroutine helper_svgrps_Nter(imol,rs,resname,sty,sty2)
!
  use iounit
  use molecule
  use polypep
  use atoms
  use sequen
  use fos
  use system
!
  implicit none
!
  integer rs,imol,j,hh,tt,ender,off,shf2,shf3,shf6,shf7
  character(3) resname
  character sty,sty2
!
  shf2 = 0
  shf3 = 0
  shf6 = 0
  shf7 = 0
  if (ua_model.gt.0) then
    shf2 = 2
    shf3 = 3
    shf6 = 6
    shf7 = 7
  end if
!
  if ((resname.eq.'ACE').OR.(resname.eq.'FOR')) then
!
!   do nothing
!
  else if (sty.eq.'P') then
!
    if (resname.eq.'UNK') then
      write(ilog,*) 'Warning. Even though unsupported residue is recognized as polypeptide, solvation &
 &groups for N-terminal residue ',rs,' are not set up (patch required).'
      return
    end if
!
!   we'll overallocate and then just use however many we need
    allocate(at(rs)%fosgrp(5))
    at(rs)%nfosgrps = 0
!
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    if ((resname.eq.'PRO').OR.(resname.eq.'HYP').OR.&
 &      (resname.eq.'PCA')) then
      if (moltermid(imol,1).eq.1) then
        hh = 3
        ender = 6
        at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PEP_PROCNT(1:3)
      else if (moltermid(imol,1).eq.2) then
        hh = 2
        ender = 5
        at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PEP_PROUNT(1:3)
      else
        write(ilog,*) 'Fatal. N-terminus type ',&
 &moltermid(imol,1),'(mol. ',imol,') is not yet supported by implici&
 &it solvent model.'
        call fexit()
      end if
    else
      if (moltermid(imol,1).eq.1) then
        hh = 4
        ender = 7
        at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PEP_CNT(1:3)
      else if (moltermid(imol,1).eq.2) then
        hh = 3
        ender = 6
        at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PEP_UNT(1:3)
      else
        write(ilog,*) 'Fatal. N-terminus type ',&
 &moltermid(imol,1),'(mol. ',imol,') is not yet supported by implici&
 &it solvent model.'
        call fexit()
      end if
    end if
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
!
    tt = 0
    do j=1,ender
      if (.NOT.((j.le.1).OR.(j.ge.5))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt)= at(rs)%bb(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt)= 1.0/(1.0*hh)
      which_fosg(at(rs)%bb(j)) = at(rs)%nfosgrps
    end do 
!
  else if (sty.eq.'N') then
!
    if (resname.eq.'UNK') then
      write(ilog,*) "Warning. Even though unsupported residue is recognized as nucleic acid, solvation &
 &groups for 5'-terminal residue ",rs," are not set up (patch required)."
      return
    end if
!
!   we'll overallocate and then just use however many we need
    allocate(at(rs)%fosgrp(8))
    at(rs)%nfosgrps = 0
!
!   first the phospho-monoester (XHPO4-): we'll just use the POO skeleton in the middle
!   for the PO4D-value and assign the difference to the terminal OH
    off = 0
    if (seqflag(rs).eq.22) then
      off = 5
      if (moltermid(imol,1).eq.1) then
        at(rs)%nfosgrps = at(rs)%nfosgrps + 1
        hh = 3
        allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
        allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
        at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
        at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_PO4D(1:3)
        tt = 0
        do j=2,5
          if (j.eq.3) cycle
          tt = tt + 1
          at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt)= at(rs)%bb(j)
          at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt)= 1.0/(1.0*hh)
          which_fosg(at(rs)%bb(j)) = at(rs)%nfosgrps
        end do
        at(rs)%nfosgrps = at(rs)%nfosgrps + 1
        hh = 2
        allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
        allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
        at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
        at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_PO4M(1:3) - FOS_NUC_PO4D(1:3)
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1)= at(rs)%bb(1)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1)= 1.0/(1.0*hh)
        which_fosg(at(rs)%bb(1)) = at(rs)%nfosgrps
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2) = at(rs)%bb(at(rs)%nbb)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2)= 1.0/(1.0*hh)
        which_fosg(at(rs)%bb(at(rs)%nbb)) = at(rs)%nfosgrps
      else
        write(ilog,*) 'Fatal. N-terminus type ',&
 &moltermid(imol,1),'(mol. ',imol,') is not yet supported by implici&
 &it solvent model.'
        call fexit()
      end if
!   for 5'-caps with a free OH that OH becomes the first solvation group
    else if (seqflag(rs).eq.24) then
      off = 1
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      hh = 2
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_RIBO1(1:3) - FOS_NUC_MTHF(1:3)
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1)= at(rs)%bb(1)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1)= 1.0/(1.0*hh)
      which_fosg(at(rs)%bb(1)) = at(rs)%nfosgrps
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2) = at(rs)%bb(at(rs)%nbb)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2)= 1.0/(1.0*hh)
      which_fosg(at(rs)%bb(at(rs)%nbb)) = at(rs)%nfosgrps
    end if
!
!   now the (rest of the) sugar: THF + decorating hydroxyls are really the underlying model
!   compounds for us
    if ((resname.eq.'DPC').OR.(resname.eq.'DPU').OR.&
 &      (resname.eq.'DPT').OR.(resname.eq.'DPA').OR.&
 &      (resname.eq.'DPG').OR.(resname.eq.'DIT').OR.&
 &      (resname.eq.'DIC').OR.(resname.eq.'DIU').OR.&
 &      (resname.eq.'DIG').OR.(resname.eq.'DIA')) then
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      hh = 13 - shf7
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_MTHF(1:3)
      tt = 0
      do j=1+off,at(rs)%nbb-1
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%bb(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%bb(j)) = at(rs)%nfosgrps
      end do
      do j=1,6-shf3
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
      end do
    else if ((resname.eq.'RPC').OR.(resname.eq.'RPU').OR.&
 &           (resname.eq.'RPT').OR.(resname.eq.'RPA').OR.&
 &           (resname.eq.'RPG').OR.(resname.eq.'RIT').OR.&
 &           (resname.eq.'RIC').OR.(resname.eq.'RIU').OR.&
 &           (resname.eq.'RIG').OR.(resname.eq.'RIA')) then
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      hh = 12 - shf6
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_MTHF(1:3)
      tt = 0
      do j=1+off,at(rs)%nbb-1
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%bb(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%bb(j)) = at(rs)%nfosgrps
      end do
      do j=1,6-shf3
        if (j.eq.4) cycle
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
      end do
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      hh = 2
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      if ((resname.eq.'RIT').OR.&
 &        (resname.eq.'RIC').OR.(resname.eq.'RIU').OR.&
 &        (resname.eq.'RIG').OR.(resname.eq.'RIA')) then
        at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_RIBO2(1:3) - FOS_NUC_RIBO1(1:3)
      else
        at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_RIBO1(1:3) - FOS_NUC_MTHF(1:3)
      end if
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1) = at(rs)%sc(4)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(4)) = at(rs)%nfosgrps
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2) = at(rs)%sc(7-shf2)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(7-shf2)) = at(rs)%nfosgrps
    else if ((resname.eq.'D5P').OR.(resname.eq.'DIB')) then
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      hh = 13 - shf7
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_MTHF(1:3)
      tt = 0
      do j=1+off,at(rs)%nbb-1
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%bb(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%bb(j)) = at(rs)%nfosgrps
      end do
      do j=1,6-shf3
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
      end do
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      hh = 2
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      if (resname.eq.'DIB') then
        at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_RIBO2(1:3) - FOS_NUC_RIBO1(1:3)
      else
        at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_RIBO1(1:3) - FOS_NUC_MTHF(1:3)
      end if
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2) = at(rs)%sc(7-shf3)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(7-shf3)) = at(rs)%nfosgrps
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1) = at(rs)%sc(8-shf3)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(8-shf3)) = at(rs)%nfosgrps
    else if ((resname.eq.'R5P').OR.(resname.eq.'RIB')) then
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      hh = 12 - shf6
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_MTHF(1:3)
      tt = 0
      do j=1+off,at(rs)%nbb-1
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%bb(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%bb(j)) = at(rs)%nfosgrps
      end do
      do j=1,6-shf3
        if (j.eq.4) cycle
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
      end do
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      hh = 2
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      if (resname.eq.'RIB') then
        at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_RIBO2(1:3) - FOS_NUC_RIBO1(1:3)
      else
        at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_RIBO1(1:3) - FOS_NUC_MTHF(1:3)
      end if
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1) = at(rs)%sc(4)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(4)) = at(rs)%nfosgrps
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2) = at(rs)%sc(7-shf2)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(7-shf2)) = at(rs)%nfosgrps
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      hh = 2
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      if (resname.eq.'RIB') then
        at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_RIBO3(1:3) - FOS_NUC_RIBO2(1:3)
      else
        at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_RIBO2(1:3) - FOS_NUC_RIBO1(1:3)
      end if
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1) = at(rs)%sc(8-shf2)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(8-shf2)) = at(rs)%nfosgrps
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2) = at(rs)%sc(9-shf2)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(9-shf2)) = at(rs)%nfosgrps
    end if
!
  else if (seqflag(rs).eq.10) then
    write(ilog,*) 'Warning. N-terminal polypeptide cap residue will &
 &have no solvation groups set up (patch required for residue ',rs,').'
    return
!
  else if (resname.eq.'UNK') then
    write(ilog,*) 'Warning. N-terminal cap residue of unrecognized polymer type will &
 &have no solvation groups set up (patch required for residue ',rs,').'
    return
!
  else
!
    write(ilog,*) 'Fatal. Implicit solvent does not support this typ&
 &e of polymer or these polypeptide caps at the moment.'
    call fexit()
!
  end if
!
end
!
!-----------------------------------------------------------------------
!
! the final routine for backbone-associated groups -> deals with residues
! that are non-terminal
!
subroutine helper_svgrps_middle(imol,rs,resname,resname2,sty,sty2)
!
  use iounit
  use molecule
  use polypep
  use atoms
  use sequen
  use fos
  use system
!
  implicit none
!
  integer rs,imol,j,hh,tt,shf2,shf3,shf6,shf7
  character(3) resname,resname2
  character sty,sty2
!
  shf2 = 0
  shf3 = 0
  shf6 = 0
  shf7 = 0
  if (ua_model.gt.0) then
    shf2 = 2
    shf3 = 3
    shf6 = 6
    shf7 = 7
  end if
!
  if (sty.eq.'P') then
!
!   we'll overallocate and then just use however many we need
    allocate(at(rs)%fosgrp(4))
    at(rs)%nfosgrps = 0
!
    if ((ci(rs-1).le.0).OR.(oi(rs-1).le.0)) then
      if (resname2.eq.'UNK') then
        write(ilog,*) 'Warning. Unusual polymer linkage due to unsupported residue leads to elimination of &
 &backbone solvation group for residue ',rs,' in helper_svgrps_middle(...).'
        at(rs)%nfosgrps = 0
        return
      else
        write(ilog,*) 'Fatal. Reference atoms missing for assumed polypeptide linkage &
 &for the backbone solvation group for residue ',rs,' in helper_svgrps_middle(...). This is a bug.'
        call fexit()
      end if
    end if
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    if ((seqflag(rs).eq.5).OR.(hni(rs).le.0)) then
      hh = 3
    else
      hh = 4
    end if
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
!
!   first add CO on the previous residue
!   note we use a shift-by-100 trick to indicate that ci(rs-1) and oi(rs-1)
!   are not part of residue rs-1's fos-group list
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1)= ci(rs-1)
    at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1)= 1.0/(1.0*hh)
    which_fosg(ci(rs-1)) = at(rs)%nfosgrps + 100
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2)= oi(rs-1)
    at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2)= 1.0/(1.0*hh)
    which_fosg(oi(rs-1)) = at(rs)%nfosgrps + 100
!   set the right ref. value
    if (seqflag(rs).eq.5) then
      if (resname2.eq.'FOR') then
        at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PEP_PROBB_FOR(1:3)
      else
        at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PEP_PROBB(1:3)
      end if
    else
      if (resname2.eq.'FOR') then
        at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) =FOS_PEP_BB_FOR(1:3)
      else
        at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PEP_BB(1:3)
      end if
    end if
!   add the NH-atoms
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(3)= ni(rs)
    at(rs)%fosgrp(at(rs)%nfosgrps)%wts(3)= 1.0/(1.0*hh)
    which_fosg(ni(rs)) = at(rs)%nfosgrps
    if ((seqflag(rs).ne.5).AND.(hni(rs).gt.0)) then
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(4)= hni(rs)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(4)= 1.0/(1.0*hh)
      which_fosg(hni(rs)) = at(rs)%nfosgrps
    end if
!
  else if (sty.eq.'N') then
!
    if (resname.eq.'UNK') then
      write(ilog,*) 'Warning. Even though unsupported residue is recognized as nucleic acid, solvation &
 &groups for residue ',rs,' are not set up (patch required).'
      return
    end if
!   we'll overallocate and then just use however many we need
    allocate(at(rs)%fosgrp(6))
    at(rs)%nfosgrps = 0
!
!   first the phospho-diester (PO4-): we'll just use the POO skeleton in the middle
!   since the estered oxygens cannot carry charge (but overlap large anyway)
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 3
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_PO4D(1:3)
    tt = 0
    do j=2,5
      if (j.eq.3) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%bb(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%bb(j)) = at(rs)%nfosgrps
    end do
!
!   now the sugar: THF + decorating hydroxyls are really the underlying model
!   compounds for us
    if ((resname.eq.'DPC').OR.(resname.eq.'DPU').OR.&
 &      (resname.eq.'DPT').OR.(resname.eq.'DPA').OR.&
 &      (resname.eq.'DPG')) then
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      hh = 13 - shf7
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_MTHF(1:3)
      tt = 0
      do j=6,at(rs)%nbb
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%bb(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%bb(j)) = at(rs)%nfosgrps
      end do
      do j=1,6-shf3
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
      end do
    else if ((resname.eq.'RPC').OR.(resname.eq.'RPU').OR.&
 &           (resname.eq.'RPT').OR.(resname.eq.'RPA').OR.&
 &           (resname.eq.'RPG')) then
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      hh = 12 - shf6
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_MTHF(1:3)
      tt = 0
      do j=6,at(rs)%nbb
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%bb(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%bb(j)) = at(rs)%nfosgrps
      end do
      do j=1,6-shf3
        if (j.eq.4) cycle
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
      end do
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      hh = 2
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_RIBO1(1:3) - FOS_NUC_MTHF(1:3)
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1) = at(rs)%sc(4)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(4)) = at(rs)%nfosgrps
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2) = at(rs)%sc(7-shf2)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(7-shf2)) = at(rs)%nfosgrps
    else if (resname.eq.'D5P') then
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      hh = 13 - shf7
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_MTHF(1:3)
      tt = 0
      do j=6,at(rs)%nbb
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%bb(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%bb(j)) = at(rs)%nfosgrps
      end do
      do j=1,6-shf3
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
      end do
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      hh = 2
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_RIBO1(1:3) - FOS_NUC_MTHF(1:3)
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2) = at(rs)%sc(7-shf3)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(7-shf3)) = at(rs)%nfosgrps
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1) = at(rs)%sc(8-shf3)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(8-shf3)) = at(rs)%nfosgrps
    else if (resname.eq.'R5P') then
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      hh = 12 - shf6
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_MTHF(1:3)
      tt = 0
      do j=6,at(rs)%nbb
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%bb(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%bb(j)) = at(rs)%nfosgrps
      end do
      do j=1,6-shf3
        if (j.eq.4) cycle
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
      end do
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      hh = 2
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_RIBO1(1:3) - FOS_NUC_MTHF(1:3)
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1) = at(rs)%sc(4)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(4)) = at(rs)%nfosgrps
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2) = at(rs)%sc(7-shf2)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(7-shf2)) = at(rs)%nfosgrps
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      hh = 2
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_RIBO2(1:3) - FOS_NUC_RIBO1(1:3)
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1) = at(rs)%sc(8-shf2)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(8-shf2)) = at(rs)%nfosgrps
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2) = at(rs)%sc(9-shf2)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(9-shf2)) = at(rs)%nfosgrps
    end if
!
  else if (resname.eq.'UNK') then
!   do nothing - report elsewhere 
!
  else
!
    write(ilog,*) 'Fatal. Implicit solvent does not support this typ&
 &e of polymer or these polypeptide caps at the moment.'
    call fexit()
!
  end if
!
end
!
!-----------------------------------------------------------------------
!
! finally, two routines for sidechain-associated groups -> first one deals 
! with peptide sidechains
!
subroutine helper_svgrps_scpep(imol,rs,resname,sty,sty2)
!
  use iounit
  use molecule
  use polypep
  use atoms
  use sequen
  use fos
  use fyoc
  use system
!
  implicit none
!
  integer rs,imol,j,hh,tt,shf,shf2,shf3,shf5,shf7,shf9,lk,shfx1,shfx4,shfx2
  character(3) resname
  character sty,sty2
!
  shf = 0
  shf2 = 0
  shf3 = 0
  shf5 = 0
  shf7 = 0
  shf9 = 0
  shfx4 = 0
  shfx2 = 0
  shfx1 = 0
  if (ua_model.gt.0) then
    shf = 1
    shf2 = 2
    shf3 = 3
    shf5 = 5
    shf7 = 7
    shf9 = 9
    if (ua_model.eq.2) then
      shfx2 = 2
      shfx4 = 4
      shfx1 = 1
    end if
  end if
!
! first deal with sc-sc-crosslinks
  if (disulf(rs).gt.0) then
    lk = crlk_idx(rs)
    if ((crosslink(lk)%itstype.eq.1).OR.(crosslink(lk)%itstype.eq.2)) then
!     disulfide linkage -> we'll assume an isotropic CH2S group
      at(rs)%nfosgrps = at(rs)%nfosgrps + 1
      hh = at(rs)%nsc - 1 + shf
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
      allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
      at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
      at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3)
      tt = 0
      do j=2-shf,at(rs)%nsc
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
      end do
    else
      write(ilog,*) 'Fatal. Encountered unsupported crosslink type in helper_svgrps_scpep(...).&
 & This is most likely an omission bug.'
      call fexit()
    end if
    return
  end if
!
  if ((resname.eq.'GLY').OR.(resname.eq.'RIB').OR.&
 &    (resname.eq.'DIB').OR.(resname.eq.'R5P').OR.&
 &    (resname.eq.'D5P').OR.(resname.eq.'UNK').OR.&
 &       ((sty.ne.'P').AND.(sty.ne.'N'))) then
!   do nothing
!
! isotropic (in a solvation sense) sidechains 
  else if ((resname.eq.'ALA').OR.(resname.eq.'VAL').OR.&
 &         (resname.eq.'NVA').OR.(resname.eq.'LEU').OR.&
 &         (resname.eq.'ILE').OR.(resname.eq.'NLE').OR.&
 &         (resname.eq.'PRO').OR.&
 &         (resname.eq.'PHE').OR.(resname.eq.'ABA')) then
!
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = at(rs)%nsc - 1 + shf
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3)
    tt = 0
    do j=2-shf,at(rs)%nsc
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!
! AIB
  else if (resname.eq.'AIB') then
!
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = at(rs)%nsc 
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3)
    tt = 0
    do j=1,at(rs)%nsc
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!
! serine and cysteine (polar sites, and almost nothing else)
  else if ((resname.eq.'CYX').OR.(resname.eq.'CYS').OR.(resname.eq.'SER')) then
!
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 2
    if ((resname.eq.'CYX').OR.((ua_model.eq.2).AND.(resname.eq.'CYS'))) hh = 1
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3)
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1)= at(rs)%sc(3-shf)
    at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1)= 1.0/(1.0*hh)
    which_fosg(at(rs)%sc(3-shf)) = at(rs)%nfosgrps
    if ((resname.ne.'CYX').AND.((ua_model.lt.2).OR.(resname.ne.'CYS'))) then
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2)= at(rs)%sc(6-shf3)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2)= 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(6-shf3)) = at(rs)%nfosgrps
    end if
!
! phosphoserine (phosphate dominates, strategy similar to 5' phosphate in nucleotides)
  else if (resname.eq.'SEP') then
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 4
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_PO4D(1:3)
    tt = 0
    do j=3-shf,at(rs)%nsc
      if (j.ge.(8-shf)) exit
      if (j.eq.(5-shf)) cycle ! terminal O
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 2
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3) - FOS_NUC_PO4D(1:3)
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1) = at(rs)%sc(5-shf)
    at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1) = 1.0/(1.0*hh)
    which_fosg(at(rs)%sc(5-shf)) = at(rs)%nfosgrps
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2) = at(rs)%sc(at(rs)%nsc)
    at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2) = 1.0/(1.0*hh)
    which_fosg(at(rs)%sc(at(rs)%nsc)) = at(rs)%nfosgrps
!
! threonine (polar site and a methyl-group, split up (the methyl group part is derived from exp. solvation
!        free energies for ethanol vs. methanol), if the sidechain is fully solvated,
!        the correct solvation free energy is recovered)
  else if (resname.eq.'THR') then
!   the OH-part
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 2
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_SER(1:3)
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1)= at(rs)%sc(3-shf)
    at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1)= 1.0/(1.0*hh)
    which_fosg(at(rs)%sc(3-shf)) = at(rs)%nfosgrps
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2)= at(rs)%sc(6-shf2)
    at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2)= 1.0/(1.0*hh)
    which_fosg(at(rs)%sc(6-shf2)) = at(rs)%nfosgrps
!   the CH3-part
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 4-shf3
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3) - FOS_SER(1:3)
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1)= at(rs)%sc(4-shf)
    at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1)= 1.0/(1.0*hh)
    which_fosg(at(rs)%sc(4-shf)) = at(rs)%nfosgrps
    tt = 1
    do j=7,at(rs)%nsc
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!
! phosphothreonine (phosphate dominates, strategy combined of THR and 5' phosphate in nucleotides)
  else if (resname.eq.'TPO') then
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 4
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_PO4D(1:3)
    tt = 0
    do j=3-shf,at(rs)%nsc
      if (j.ge.(9-shf)) exit
      if ((j.eq.(4-shf)).OR.(j.eq.(6-shf))) cycle ! CG, terminal O
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!   terminal OH
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 2
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_SEP(1:3) - FOS_NUC_PO4D(1:3)
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1) = at(rs)%sc(6-shf)
    at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1) = 1.0/(1.0*hh)
    which_fosg(at(rs)%sc(6-shf)) = at(rs)%nfosgrps
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2) = at(rs)%sc(at(rs)%nsc)
    at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2) = 1.0/(1.0*hh)
    which_fosg(at(rs)%sc(at(rs)%nsc)) = at(rs)%nfosgrps
!   the CH3-part
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 4-shf3
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3) - FOS_SEP(1:3)
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1)= at(rs)%sc(4-shf)
    at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1)= 1.0/(1.0*hh)
    which_fosg(at(rs)%sc(4-shf)) = at(rs)%nfosgrps
    tt = 1
    if (ua_model.eq.0) then
      do j=10,12
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
      end do
    end if
!
! tyrosine (polar site and toluene-system, same strategy as for threonine)
  else if (resname.eq.'TYR') then
!   the OH-part
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 2
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3) - FOS_PHE(1:3)
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1)= at(rs)%sc(9-shf)
    at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1)= 1.0/(1.0*hh)
    which_fosg(at(rs)%sc(9-shf)) = at(rs)%nfosgrps
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2)= at(rs)%sc(16-shf3-shfx4)
    at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2)= 1.0/(1.0*hh)
    which_fosg(at(rs)%sc(16-shf3-shfx4)) = at(rs)%nfosgrps
!   the toluene part
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = at(rs)%nsc - 3 + shf
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PHE(1:3)
    tt = 0
    do j=2-shf,at(rs)%nsc
      if ((j.eq.(9-shf)).OR.(j.eq.(16-shf3-shfx4))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!
! tyrosine(-) (charge site and toluene-system, same strategy as for tyrosine)
  else if (resname.eq.'TYO') then
!   the O(-)-part
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 1
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3) - FOS_PHE(1:3)
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1)= at(rs)%sc(9-shf)
    at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1)= 1.0/(1.0*hh)
    which_fosg(at(rs)%sc(9-shf)) = at(rs)%nfosgrps
!   the toluene part
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = at(rs)%nsc - 2 + shf
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PHE(1:3)
    tt = 0
    do j=2-shf,at(rs)%nsc
      if (j.eq.(9-shf)) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!
! phosphotyrosine (polar site and toluene-system, same strategy as for tyrosine)
  else if (resname.eq.'PTR') then
!   phosphate
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 4
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NUC_PO4D(1:3)
    tt = 0
    do j=9-shf,at(rs)%nsc
      if (j.ge.(14-shf)) exit
      if (j.eq.(11-shf)) cycle ! terminal O
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!   the OH-part
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 2
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3) - FOS_PHE(1:3) - FOS_NUC_PO4D(1:3)
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1)= at(rs)%sc(11-shf)
    at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1)= 1.0/(1.0*hh)
    which_fosg(at(rs)%sc(11-shf)) = at(rs)%nfosgrps
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2)= at(rs)%sc(at(rs)%nsc)
    at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2)= 1.0/(1.0*hh)
    which_fosg(at(rs)%sc(at(rs)%nsc)) = at(rs)%nfosgrps
!   the toluene part
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = at(rs)%nsc - 7 + shf
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PHE(1:3)
    tt = 0
    do j=2-shf,at(rs)%nsc-1
      if ((j.ge.(9-shf)).AND.(j.le.(13-shf))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!
! methionine: a difficult case, the methyl-ethyl thioether has a (slightly) favorable solvation FE
!         while the underlying butane sidechain has a slightly positive SFE; we will assign
!         the difference to the sulphur atom and the butane value to the rest
   else if (resname.eq.'MET') then
!   the sulphur atom
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 1
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3) - FOS_NLE(1:3)
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1)= at(rs)%sc(4-shf)
    at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1)= 1.0/(1.0*hh)
    which_fosg(at(rs)%sc(4-shf)) = at(rs)%nfosgrps
!   the rest
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = at(rs)%nsc - 2 + shf
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NLE(1:3)
    tt = 0
    do j=2-shf,at(rs)%nsc
      if ((j.eq.(4-shf))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!
! tryptophane: another difficult case, the indol ring has a favorable SFE, but the decompisition into
!          the polar (but not excessively so) NH-part and the large ring system is not straightforward
!          we will use the SFE of naphthalene to get a value for the NH-free indol-like system (makes a ~60/40 split)
  else if (resname.eq.'TRP') then
!   the polar N-H
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 2
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3) - FOS_NAP(1:3)
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1)= at(rs)%sc(6-shf)
    at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1)= 1.0/(1.0*hh)
    which_fosg(at(rs)%sc(6-shf)) = at(rs)%nfosgrps
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2)= at(rs)%sc(15-shf3-shfx1)
    at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2)= 1.0/(1.0*hh)
    which_fosg(at(rs)%sc(15-shf3-shfx1)) = at(rs)%nfosgrps
!   the rest
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = at(rs)%nsc - 3 + shf
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_NAP(1:3)
    tt = 0
    do j=2-shf,at(rs)%nsc
      if ((j.eq.(6-shf)).OR.(j.eq.(15-shf3-shfx1))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!
! histidine N_epsilon-protonated: another difficult case, the ring is polar overall but
!                             the primary hydration sites are clearly the two
!                             nitrogens; we will use the four atoms NH-C-N for now
!                             unlike TRP/TYR, there's no clear polar domain to single out
!                             and no base compound to compare with 
  else if (resname.eq.'HIE') then
!   the polar N-C-NH
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 4
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3)
    tt = 0
    do j=2-shf,at(rs)%nsc
      if (.NOT.((j.eq.(4-shf)).OR.(j.eq.(6-shf)).OR.&
 &              (j.eq.(7-shf)).OR.(j.eq.(12-shf3-shfx2))))  cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!
! histidine N_delta-protonated: see above
  else if (resname.eq.'HID') then
!   the polar NH-C-N
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 4
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3)
    tt = 0
    do j=2-shf,at(rs)%nsc
      if (.NOT.((j.eq.(4-shf)).OR.(j.eq.(6-shf)).OR.&
 &              (j.eq.(7-shf)).OR.(j.eq.(10-shf3)))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!
! histidine doubly protonated: see above
  else if (resname.eq.'HIP') then
!   the polar NH-C-NH
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 5
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3)
    tt = 0
    do j=2-shf,at(rs)%nsc
      if (.NOT.((j.eq.(4-shf)).OR.(j.eq.(6-shf)).OR.&
 &              (j.eq.(7-shf)).OR.(j.eq.(10-shf3).OR.(j.eq.(13-shf3-shfx2))))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!
! asparagine (let's focus attention on the polar sites and ignore the single methylene group)
  else if (resname.eq.'ASN') then
!
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 5
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3)
    tt = 0
    do j=3-shf,at(rs)%nsc
      if ((ua_model.eq.0).AND.((j.eq.6).OR.(j.eq.7))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!
! (protonated) aspartic acid (same strategy as ASN)
  else if (resname.eq.'ASH') then
!
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 4
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3)
    tt = 0
    do j=3-shf,at(rs)%nsc
      if ((ua_model.eq.0).AND.((j.eq.6).OR.(j.eq.7))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!
! glutamine: the difference between acetamide and propionamide is ~0.3kcal/mol, let's assign that
!            small difference to the C-beta methylene
  else if (resname.eq.'GLN') then
!   the amide
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 5
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_ASN(1:3)
    tt = 0
    do j=2-shf,at(rs)%nsc
      if (.NOT.((j.eq.(4-shf)).OR.(j.eq.(5-shf)).OR.(j.eq.(6-shf)).OR.&
 &                             (j.eq.(11-shf5)).OR.(j.eq.(12-shf5)))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!   the Cbeta-methylene
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 3-shf2
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3) - FOS_ASN(1:3)
    tt = 0
    do j=2-shf,at(rs)%nsc
      if (.NOT.((j.eq.(2-shf)).OR.(j.eq.7).OR.(j.eq.8))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
      if ((ua_model.gt.0).AND.(j.ge.1)) exit
    end do

! (protonated) glutamic acid: the difference between acetic acid and propionic acid is ~0.2kcal/mol,
!                             let's assign that small difference to the C-beta methylene
  else if (resname.eq.'GLH') then
!   the acid
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 4
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_ASH(1:3)
    tt = 0
    do j=2-shf,at(rs)%nsc
      if (.NOT.((j.eq.(4-shf)).OR.(j.eq.(5-shf)).OR.(j.eq.(6-shf)).OR.&
 &                             (j.eq.(11-shf5)))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!   the Cbeta-methylene
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 3-shf2
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3) - FOS_ASH(1:3)
    tt = 0
    do j=2-shf,at(rs)%nsc
      if (.NOT.((j.eq.(2-shf)).OR.(j.eq.7).OR.(j.eq.8))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
      if ((ua_model.gt.0).AND.(j.ge.1)) exit
    end do
!
! aspartate (charge localized to COO for solvation purposes, overwhelming -> rest neglected)
  else if (resname.eq.'ASP') then
!   the carboxyl group
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 3
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3)
    tt = 0
    do j=2-shf,at(rs)%nsc
      if (.NOT.((j.eq.(3-shf)).OR.(j.eq.(4-shf)).OR.(j.eq.(5-shf)))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!
! glutamate (same strategy as aspartate)
  else if (resname.eq.'GLU') then
!   the carboxyl group
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 3
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3)
    tt = 0
    do j=2-shf,at(rs)%nsc
      if (.NOT.((j.eq.(6-shf)).OR.(j.eq.(4-shf)).OR.(j.eq.(5-shf)))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!
! lysine: we will use the value for the full analog localized on the amino group
  else if (resname.eq.'LYS') then
!   the ammonium ion
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 4
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3)
    tt = 0
    do j=2-shf,at(rs)%nsc
      if (.NOT.((j.eq.(6-shf)).OR.(j.ge.(at(rs)%nsc-2)))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!
! deprotonated lysine: difference between full value and propane to amine group and propane to CB-CD
  else if (resname.eq.'LYD') then
!   the amine group
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 3
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3) - FOS_PRP
    tt = 0
    do j=2-shf,at(rs)%nsc
      if (.NOT.((j.eq.(6-shf)).OR.(j.ge.(at(rs)%nsc-1)))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 9
    if (ua_model.gt.0) hh = 3 
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PRP
    tt = 0
    do j=2-shf,12-shf9
      if ((j.eq.(6-shf)).OR.(j.eq.(5-shf)).OR.(j.ge.(at(rs)%nsc-1))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!
! acetylated lysine: difference between full value and propane to amide (NMA) group and propane to C3 part
  else if (resname.eq.'KAC') then
!   the amide group: as usual focus on the polar atoms
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 4
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3) - FOS_PRP
    at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1:hh) = 1.0/(1.0*hh)
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1) = at(rs)%sc(6-shf)
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2) = at(rs)%sc(7-shf)
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(3) = at(rs)%sc(9-shf) 
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(4) = at(rs)%sc(18-shf9)
    which_fosg(at(rs)%sc(6-shf)) = at(rs)%nfosgrps
    which_fosg(at(rs)%sc(7-shf)) = at(rs)%nfosgrps
    which_fosg(at(rs)%sc(9-shf)) = at(rs)%nfosgrps
    which_fosg(at(rs)%sc(18-shf9)) = at(rs)%nfosgrps
!   the aliphatic tail (CB-CD)
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 9
    if (ua_model.gt.0) hh = 3 
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PRP
    tt = 0
    do j=2-shf,min(at(rs)%nsc,15-shf)
      if (.NOT.((j.le.(4-shf)).OR.(j.ge.10))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!
! monomethylated lysine: difference between full value and propane to ammonium ion and propane to CB-CD
  else if (resname.eq.'KM1') then
!   the ammonium ion (just carbon from methyl)
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 4
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3) - FOS_PRP
    at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1:hh) = 1.0/(1.0*hh)
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1) = at(rs)%sc(6-shf)
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2) = at(rs)%sc(7-shf)
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(3) = at(rs)%sc(16-shf9) 
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(4) = at(rs)%sc(17-shf9)
    which_fosg(at(rs)%sc(6-shf)) = at(rs)%nfosgrps
    which_fosg(at(rs)%sc(7-shf)) = at(rs)%nfosgrps
    which_fosg(at(rs)%sc(16-shf9)) = at(rs)%nfosgrps
    which_fosg(at(rs)%sc(17-shf9)) = at(rs)%nfosgrps
!   the aliphatic tail (CB-CD)
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 9
    if (ua_model.gt.0) hh = 3 
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PRP
    tt = 0
    do j=2-shf,min(at(rs)%nsc-2,13-shf)
      if (.NOT.((j.le.(4-shf)).OR.(j.ge.8))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!
! dimethylated lysine: difference between full value and propane to ammonium ion and propane to CB-CD
  else if (resname.eq.'KM2') then
!   the ammonium ion (just carbons from methyls)
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 4
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3) - FOS_PRP
    at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1:hh) = 1.0/(1.0*hh)
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1) = at(rs)%sc(6-shf)
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2) = at(rs)%sc(7-shf)
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(3) = at(rs)%sc(8-shf) 
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(4) = at(rs)%sc(17-shf9)
    which_fosg(at(rs)%sc(6-shf)) = at(rs)%nfosgrps
    which_fosg(at(rs)%sc(7-shf)) = at(rs)%nfosgrps
    which_fosg(at(rs)%sc(8-shf)) = at(rs)%nfosgrps
    which_fosg(at(rs)%sc(17-shf9)) = at(rs)%nfosgrps
!   the aliphatic tail (CB-CD)
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 9
    if (ua_model.gt.0) hh = 3 
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PRP
    tt = 0
    do j=2-shf,min(at(rs)%nsc-1,14-shf)
      if (.NOT.((j.le.(4-shf)).OR.(j.ge.9))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!
! trimethylated lysine: simply assume the entire TMA model compound; ignore rest
  else if (resname.eq.'KM3') then
!   the ammonium ion (just carbons from methyls)
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 5
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3)
    at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1:hh) = 1.0/(1.0*hh)
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1) = at(rs)%sc(5-shf)
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2) = at(rs)%sc(6-shf)
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(3) = at(rs)%sc(7-shf)
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(4) = at(rs)%sc(8-shf)
    at(rs)%fosgrp(at(rs)%nfosgrps)%ats(5) = at(rs)%sc(9-shf)
    which_fosg(at(rs)%sc(5-shf)) = at(rs)%nfosgrps
    which_fosg(at(rs)%sc(6-shf)) = at(rs)%nfosgrps
    which_fosg(at(rs)%sc(7-shf)) = at(rs)%nfosgrps
    which_fosg(at(rs)%sc(8-shf)) = at(rs)%nfosgrps
    which_fosg(at(rs)%sc(9-shf)) = at(rs)%nfosgrps
!
! ornithine: we will use the value for the full analog localized on the amino group
  else if (resname.eq.'ORN') then
!   the ammonium ion
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 4
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3)
    tt = 0
    do j=2-shf,at(rs)%nsc
      if (.NOT.((j.eq.(5-shf)).OR.(j.ge.(at(rs)%nsc-2)))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!
! DAB: we will use the value for the full analog localized on the amino group
  else if (resname.eq.'DAB') then
!   the ammonium ion
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 4
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3)
    tt = 0
    do j=2-shf,at(rs)%nsc
      if (.NOT.((j.eq.(4-shf)).OR.(j.ge.(at(rs)%nsc-2)))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!
! arginine: the solvation is overwhelmingly driven by the delocalized charge in the guanidino
!       group, for the sake of simplicity the aliphatic connector (3C) is ignored
  else if (resname.eq.'ARG') then
!   the five polar hydrogens (15-19) and the four heavy atoms of the guano-group (5-8)
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 9
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3)
    tt = 0
    do j=2-shf,at(rs)%nsc
      if (.NOT.((j.ge.15-shf7).OR.((j.ge.5-shf).AND.(j.le.8-shf)))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!
  else if (sty.eq.'N') then
!
!    do nothing yet (see below!!!)
!
! all other residues are not yet supported (or never will be)
  else
    write(ilog,*) 'Fatal. Implicit solvent does not support sidechai&
 &n of residue ',resname,' at the moment. Check back later.'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
! second one deals with nucleic acid sidechains
!
subroutine helper_svgrps_scnuc(imol,rs,resname,sty,sty2)
!
  use iounit
  use molecule
  use polypep
  use atoms
  use sequen
  use fos
  use system
!
  implicit none
!
  integer rs,imol,j,hh,tt,off,shf3,shfx1,shfx2
  character(3) resname
  character sty,sty2
!
  shf3 = 0
  shfx1 = 0
  shfx2 = 0
  if (ua_model.gt.0) then
    shf3 = 3
    if (ua_model.eq.2) then
      shfx2 = 2
      shfx1 = 1
    end if
  end if
!
!
! now on to nucleic acid bases:
! we have a very fundamental problem here in that the glycosyl-bond eliminates a H-donor on both
! ends (the sugar AND the base) -> this is not amenable to the same set of assumptions we have
! with parsing polypeptides into model compounds where the sidechain link is always aliphatic and
! the backbone link is assumed to be through the alpha-carbon -> also aliphatic
! the corresponding model compounds are therefore not really ribose and pyr/pur, but instead should
! probably be various forms of hydroxy-THFs and N->C substitutes of pyr/pur
! to be clear, though, the problem is in the parameters, not the implementation
!
! uracil is very favorably solvated (ab initio calcs), but in the glycolysated form loses on of the
! two H-bond donors (which can have a huge effect as told by the difference A->G). we'll assign
! the favorable part minus benzene to the NCONHCO and benzene to the CHCH
  if ((resname.eq.'DIU').OR.(resname.eq.'RIU').OR.&
 &         (resname.eq.'DPU').OR.(resname.eq.'RPU')) then
    off = 0
    if ((resname.eq.'RPU').OR.(resname.eq.'RIU')) off = 1
!   the polar part
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    if (ua_model.gt.0) then
      off = -3
      if ((resname.eq.'RPU').OR.(resname.eq.'RIU')) off = -1
    end if
    hh = 7
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = freesolv(rs,1:3) - FOS_BEN(1:3)
    tt = 0
    do j=7+off,at(rs)%nsc
      if ((j.eq.(9+off)).OR.(j.eq.(12+off))) cycle
      if ((ua_model.lt.2).AND.(j.eq.(15+off)).OR.(j.eq.(17+off))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!   the CHCH-part
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 4
    if (ua_model.eq.2) hh = 2
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_BEN(1:3)
    tt = 0
    if (ua_model.lt.2) then
      do j=9+off,at(rs)%nsc
        if (.NOT.((j.eq.(9+off)).OR.(j.eq.(12+off)).OR.&
 &         (j.eq.(15+off)).OR.(j.eq.(17+off)))) cycle
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
      end do
    else
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1) = at(rs)%sc(9+off)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(9+off)) = at(rs)%nfosgrps
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2) = at(rs)%sc(12+off)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(12+off)) = at(rs)%nfosgrps
    end if
!
! thymine is less favorably solvated (ab initio calcs), the difference is
! assigned to the methyl group
  else if ((resname.eq.'DIT').OR.(resname.eq.'RIT').OR.&
 &         (resname.eq.'DPT').OR.(resname.eq.'RPT')) then
    off = 0
    if ((resname.eq.'RPT').OR.(resname.eq.'RIT')) off = 1
    if (ua_model.gt.0) then
      off = -3
      if ((resname.eq.'RPT').OR.(resname.eq.'RIT')) off = -1
    end if
!   the polar part
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 7
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_XPU(1:3) - FOS_BEN(1:3)
    tt = 0
    do j=7+off,at(rs)%nsc-3+shf3
      if ((j.eq.(9+off)).OR.(j.eq.(12+off)).OR.(j.eq.(14+off))) cycle
      if ((ua_model.lt.2).AND.(j.eq.(16+off))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!   the aromatic CH-part
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 3
    if (ua_model.eq.2) hh = 2
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_BEN(1:3)
    tt = 0
    if (ua_model.lt.2) then
      do j=9+off,at(rs)%nsc-3+shf3
        if (.NOT.((j.eq.(9+off)).OR.(j.eq.(12+off)).OR.&
 &         (j.eq.(16+off)))) cycle
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
      end do
    else
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(1) = at(rs)%sc(9+off)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(1) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(9+off)) = at(rs)%nfosgrps
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(2) = at(rs)%sc(12+off)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(2) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(12+off)) = at(rs)%nfosgrps
    end if
!   the methyl group
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 4 - shf3
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_XPT(1:3) - FOS_XPU(1:3)
    tt = 0
    do j=14+off,at(rs)%nsc
      if (.NOT.((j.eq.(14+off)).OR.(j.gt.(at(rs)%nsc-3)))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
      if ((ua_model.gt.0).AND.(j.ge.(14+off))) exit
    end do
!
! cytosine is eminently more favorable than uracil: we'll assign the extra to the free aniline-style amine 
  else if ((resname.eq.'DIC').OR.(resname.eq.'RIC').OR.&
 &         (resname.eq.'DPC').OR.(resname.eq.'RPC')) then
    off = 0
    if ((resname.eq.'RPC').OR.(resname.eq.'RIC')) off = 1
    if (ua_model.gt.0) then
      off = -3
      if ((resname.eq.'RPC').OR.(resname.eq.'RIC')) off = -1
    end if
!   the polar part (only NCON this time)
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 4
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_XPU(1:3) - FOS_BEN(1:3)
    tt = 0
    do j=7+off,11+off
      if (j.eq.(9+off)) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!   the aromatic CH-part (CCHCH this time)
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 5
    if (ua_model.eq.2) hh = 3
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_BEN(1:3)
    tt = 0
    if (ua_model.lt.2) then
      do j=9+off,at(rs)%nsc
        if (.NOT.((j.eq.(9+off)).OR.(j.eq.(12+off)).OR.(j.eq.(13+off)).OR.&
 &              ((j.eq.(15+off))).OR.(j.eq.(16+off)))) cycle
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
      end do
    else
      do j=9+off,at(rs)%nsc
        if (.NOT.((j.eq.(9+off)).OR.(j.eq.(12+off)).OR.(j.eq.(13+off)))) cycle
        tt = tt + 1
        at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
        at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
        which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
      end do
    end if
!   the amine group
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 3
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_XPC(1:3) - FOS_XPU(1:3)
    tt = 0
    do j=14+off,at(rs)%nsc
      if (.NOT.((j.eq.(14+off)).OR.(j.gt.(at(rs)%nsc-2)))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!
! adenine has no oxygen-based acceptor, but a very good donor. the net
! FOS is on par with that of uracil/thymine, but significantly less than
! guanine or cytosine  
!
  else if ((resname.eq.'DIA').OR.(resname.eq.'RIA').OR.&
 &         (resname.eq.'DPA').OR.(resname.eq.'RPA')) then
    off = 0
    if ((resname.eq.'RPA').OR.(resname.eq.'RIA')) off = 1
    if (ua_model.gt.0) then
      off = -3
      if ((resname.eq.'RPA').OR.(resname.eq.'RIA')) off = -1
    end if
!   the polar part
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 11 - shfx2
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PUR(1:3)
    tt = 0
    do j=7+off,at(rs)%nsc-2
      if (j.eq.(16+off)) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!   the NH2-group
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 3
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_XPA(1:3) - FOS_PUR(1:3)
    tt = 0
    do j=at(rs)%nsc-4+shfx2,at(rs)%nsc
      if (ua_model.lt.2) then
        if ((j.eq.(at(rs)%nsc-3)).OR.(j.eq.(at(rs)%nsc-2))) cycle
      end if
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!
! finally, guanine has both an oxygen-based acceptor and an extra donor. with respect
! to ADE we'll assign the difference to the "new" CONH stretch  
!
  else if ((resname.eq.'DIG').OR.(resname.eq.'RIG').OR.&
 &         (resname.eq.'DPG').OR.(resname.eq.'RPG')) then
    off = 0
    if ((resname.eq.'RPG').OR.(resname.eq.'RIG')) off = 1
    if (ua_model.gt.0) then
      off = -3
      if ((resname.eq.'RPG').OR.(resname.eq.'RIG')) off = -1
    end if
!   the polar part
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 8 - shfx1
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_PUR(1:3)
    tt = 0
    do j=7+off,at(rs)%nsc-3
      if ((j.eq.(14+off)).OR.(j.eq.(15+off)).OR.&
 &        (j.eq.(16+off)).OR.(j.eq.(17+off))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!   the NH2-group
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 3
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_XPA(1:3) - FOS_PUR(1:3)
    tt = 0
    do j=at(rs)%nsc-5+shfx1,at(rs)%nsc
      if ((j.eq.(at(rs)%nsc-4+shfx1)).OR.(j.eq.(at(rs)%nsc-2))) cycle
      if ((ua_model.lt.2).AND.(j.eq.(at(rs)%nsc-3))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
!   the CONH-unit
    at(rs)%nfosgrps = at(rs)%nfosgrps + 1
    hh = 4
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%ats(hh))
    allocate(at(rs)%fosgrp(at(rs)%nfosgrps)%wts(hh))
    at(rs)%fosgrp(at(rs)%nfosgrps)%nats = hh
    at(rs)%fosgrp(at(rs)%nfosgrps)%val(1:3) = FOS_XPG(1:3) - FOS_XPA(1:3)
    tt = 0
    do j=14+off,at(rs)%nsc-2
      if (j.eq.(16+off)) cycle
      if ((ua_model.lt.2).AND.(j.eq.(18+off))) cycle
      tt = tt + 1
      at(rs)%fosgrp(at(rs)%nfosgrps)%ats(tt) = at(rs)%sc(j)
      at(rs)%fosgrp(at(rs)%nfosgrps)%wts(tt) = 1.0/(1.0*hh)
      which_fosg(at(rs)%sc(j)) = at(rs)%nfosgrps
    end do
  else if ((resname.eq.'RIB').OR.&
 &         (resname.eq.'DIB').OR.(resname.eq.'R5P').OR.&
 &         (resname.eq.'D5P').OR.&
 &         ((sty.ne.'P').AND.(sty.ne.'N')).OR.&
 &         (sty.eq.'P')) then
!
!    do nothing (see above!!!!)
!
! all other residues are not yet supported (or never will be)
  else
    write(ilog,*) 'Fatal. Implicit solvent does not support sidechai&
 &n of residue ',resname,' at the moment. Check back later.'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
! this helper routine transfer patch groups into real groups
!
subroutine helper_svgrps_patch(nps,nogps,tmpmat,tmpmat2,tmpmat3)
!
  use iounit
  use molecule
  use polypep
  use atoms
  use sequen
  use fos
  use fyoc
  use system
!
  implicit none
!
  integer tmpmat(n,2),tmpmat3(n,2),nps,nogps,i,j,k,imxfg,fs,rrr,rs,jj,kk
  RTYPE tmpmat2(n,4)
  integer, ALLOCATABLE:: grpsprs(:),ptchprs(:),restlst(:)
  type(t_fosgrp), ALLOCATABLE:: fosbu(:,:)
  logical allmtch
!
  i = 1
  allocate(restlst(n))
!
  do while (i.le.nps)
    k = tmpmat(i,2)
    j = 0 
    rrr = tmpmat(i,1)
    imxfg = i
    do while (tmpmat(i,2).eq.k)
      j = j + 1
      restlst(j) = tmpmat(i,1)
      i = i + 1
      if (i.gt.nps) exit
    end do
    allmtch = .true.
    fs = which_fosg(rrr)
    if (fs.le.0) cycle
    rs = atmres(rrr)
    if (fs.gt.100) then
      fs = fs - 100
      rs = rs + 1
    end if
    do k=1,j
      if (tmpmat2(imxfg+k-1,1).le.0.0) then
        allmtch = .false.
        exit
      end if
      if (which_fosg(restlst(k)).gt.100) then
        if ((which_fosg(restlst(k))-100).ne.fs) then
          allmtch = .false.
          exit
        end if
      else
        if (which_fosg(restlst(k)).ne.fs) then
          allmtch = .false.
          exit
        end if
      end if
    end do
    if ((allmtch.EQV..true.).AND.((i-imxfg).eq.at(rs)%fosgrp(fs)%nats)) then
!     this is the case of only atomic weights and/or rFOS changing -> change and remove from lists
      at(rs)%fosgrp(fs)%val(1:3) = tmpmat2(imxfg,2:4)
      do k=1,at(rs)%fosgrp(fs)%nats
        do rrr=1,j
          if (at(rs)%fosgrp(fs)%ats(k).eq.restlst(rrr)) then
            at(rs)%fosgrp(fs)%wts(k) = tmpmat2(rrr+imxfg-1,1)
            exit
          end if
        end do
      end do
      tmpmat(imxfg:(i-1),1) = 0 ! eliminate atom list for processed group
      do k=1,nogps
        if ((tmpmat3(k,1).eq.rs).AND.(tmpmat3(k,2).eq.fs)) then
          tmpmat3(k,1) = -tmpmat3(k,1) ! flag old group as processed
          exit
        end if
      end do
    end if
  end do
  k = 0
! shorten lists if possible
  do i=1,nps
    if (tmpmat(i,1).gt.0) then
      k = k + 1
      tmpmat(k,:) = tmpmat(i,:)
      tmpmat2(k,:) = tmpmat2(i,:)
    end if
  end do
  nps = k
!
  allocate(grpsprs(nseq))
  allocate(ptchprs(nseq))
  imxfg = maxval(at(1:nseq)%nfosgrps) + 1
  allocate(fosbu(nseq,imxfg))
!
  grpsprs(:) = at(1:nseq)%nfosgrps
  do i=1,nogps
    if (tmpmat3(i,1).lt.0) cycle
    if (allocated(fosbu(tmpmat3(i,1),1)%ats).EQV..false.) then
      do j=1,at(tmpmat3(i,1))%nfosgrps
        fosbu(tmpmat3(i,1),j)%nats = at(tmpmat3(i,1))%fosgrp(j)%nats
        fosbu(tmpmat3(i,1),j)%val(:) = at(tmpmat3(i,1))%fosgrp(j)%val(:)
        allocate(fosbu(tmpmat3(i,1),j)%wts(at(tmpmat3(i,1))%fosgrp(j)%nats))
        allocate(fosbu(tmpmat3(i,1),j)%ats(at(tmpmat3(i,1))%fosgrp(j)%nats))
        fosbu(tmpmat3(i,1),j)%wts(:) = at(tmpmat3(i,1))%fosgrp(j)%wts(:)
        fosbu(tmpmat3(i,1),j)%ats(:) = at(tmpmat3(i,1))%fosgrp(j)%ats(:)
      end do
    end if
    if (tmpmat3(i,2).eq.0) cycle
    which_fosg(at(tmpmat3(i,1))%fosgrp(tmpmat3(i,2))%ats(1:at(tmpmat3(i,1))%fosgrp(tmpmat3(i,2))%nats)) = 0
    grpsprs(tmpmat3(i,1)) = grpsprs(tmpmat3(i,1)) - 1
    at(tmpmat3(i,1))%fosgrp(tmpmat3(i,2))%nats = 0
    deallocate(at(tmpmat3(i,1))%fosgrp(tmpmat3(i,2))%ats)
    deallocate(at(tmpmat3(i,1))%fosgrp(tmpmat3(i,2))%wts)
  end do
!
  ptchprs(:) = 0
  i = 1
 22 format(a,i6,a,g10.5,a)
  do while (i.le.nps)
    j = tmpmat(i,2)
    k = atmres(tmpmat(i,1))
    do while (tmpmat(i,2).eq.j)
      rrr = min(atmres(tmpmat(i,1)),k)
      fs = max(atmres(tmpmat(i,1)),k)
      i = i + 1
      if (i.gt.nps) exit
    end do
    if (abs(rrr-fs).gt.1) then
      write(ilog,22) 'Fatal. Atoms constituting patched FOS group labeled ',j,' span more than two residues &
 &(first is ',rrr,'; last is ',fs,'). This is currently not supported. Please adjust input file.'
      call fexit()     
    end if
    if (molofrs(rrr).ne.molofrs(fs)) then
      write(ilog,22) 'Fatal. Atoms constituting patched FOS group labeled ',j,' span two molecules &
 &(first is ',molofrs(rrr),'; last is ',molofrs(fs),'). This is currently not supported. Please adjust input file.'
      call fexit()     
    end if
    ptchprs(fs) = ptchprs(fs) + 1
  end do
!
  i = 1
  do while (i.le.nps)
    k = tmpmat(i,2)
    j = atmres(tmpmat(i,1))
    rrr = i
    do while (tmpmat(rrr,2).eq.k)
      if (atmres(tmpmat(rrr,1)).gt.j) j = atmres(tmpmat(rrr,1))
      rrr = rrr + 1
      if (rrr.gt.nps) exit
    end do
!   this control is used to make sure this only done once for each affected residue
!   here, we recreate the preserved FOS groups from the backup array
    if (ptchprs(j).gt.0) then
      rrr = 0
      do fs=1,at(j)%nfosgrps
        if (at(j)%fosgrp(fs)%nats.gt.0) then
          rrr = rrr + 1
          restlst(rrr) = fs
        end if
      end do
      kk = 0
      do fs=1,at(j)%nfosgrps
        if (allocated(at(j)%fosgrp(fs)%ats).EQV..true.) then
          kk = kk + 1
          do jj=1,nogps
            if ((fs.eq.tmpmat3(jj,2)).AND.(tmpmat3(jj,1).eq.-j)) then
              tmpmat3(jj,2) = kk
            end if
          end do
          which_fosg(at(j)%fosgrp(fs)%ats(1:at(j)%fosgrp(fs)%nats)) = 0
          deallocate(at(j)%fosgrp(fs)%ats)
        end if
        if (allocated(at(j)%fosgrp(fs)%wts).EQV..true.) deallocate(at(j)%fosgrp(fs)%wts)
        at(j)%fosgrp(fs)%nats = 0
      end do
      if (allocated(at(j)%fosgrp).EQV..true.) deallocate(at(j)%fosgrp)
      allocate(at(j)%fosgrp(rrr+ptchprs(j)))
      at(j)%nfosgrps = rrr
      do fs=1,rrr
        at(j)%fosgrp(fs)%val(:) = fosbu(j,restlst(fs))%val(:)
        at(j)%fosgrp(fs)%nats = fosbu(j,restlst(fs))%nats
        allocate(at(j)%fosgrp(fs)%ats(at(j)%fosgrp(fs)%nats))
        allocate(at(j)%fosgrp(fs)%wts(at(j)%fosgrp(fs)%nats))
        at(j)%fosgrp(fs)%ats(:) = fosbu(j,restlst(fs))%ats(:)
        at(j)%fosgrp(fs)%wts(:) = fosbu(j,restlst(fs))%wts(:)
        do jj=1,at(j)%fosgrp(fs)%nats
          if (atmres(at(j)%fosgrp(fs)%ats(jj)).eq.j) then
            which_fosg(at(j)%fosgrp(fs)%ats(jj)) = fs
          else if (atmres(at(j)%fosgrp(fs)%ats(jj)).eq.(j-1)) then
            which_fosg(at(j)%fosgrp(fs)%ats(jj)) = fs+100
          else
            write(ilog,*) 'Fatal. This is a bug in helper_svgrps_patch(...). Please report.'
            call fexit()
          end if
        end do
      end do
      ptchprs(j) = 0
    end if
    rrr = i
!   now add the new group
    at(j)%nfosgrps = at(j)%nfosgrps + 1
    at(j)%fosgrp(at(j)%nfosgrps)%nats = 0
    at(j)%fosgrp(at(j)%nfosgrps)%val(1:3) = tmpmat2(i,2:4)
    do while (tmpmat(rrr,2).eq.k)
      if (tmpmat2(rrr,1).gt.0.0) then
        at(j)%fosgrp(at(j)%nfosgrps)%nats = at(j)%fosgrp(at(j)%nfosgrps)%nats + 1
      end if
      rrr = rrr + 1
      if (rrr.gt.nps) exit
    end do
    allocate(at(j)%fosgrp(at(j)%nfosgrps)%ats(at(j)%fosgrp(at(j)%nfosgrps)%nats))
    allocate(at(j)%fosgrp(at(j)%nfosgrps)%wts(at(j)%fosgrp(at(j)%nfosgrps)%nats))
    rrr = i
    fs = 0
    do while (tmpmat(i,2).eq.k)
      if (tmpmat2(i,1).gt.0.0) then
        fs = fs + 1
        at(j)%fosgrp(at(j)%nfosgrps)%ats(fs) = tmpmat(i,1)
        at(j)%fosgrp(at(j)%nfosgrps)%wts(fs) = tmpmat2(i,1)
      end if
      i = i + 1
      if (i.gt.nps) exit
    end do
    do jj=1,at(j)%fosgrp(at(j)%nfosgrps)%nats
      if (atmres(at(j)%fosgrp(at(j)%nfosgrps)%ats(jj)).eq.j) then
        which_fosg(at(j)%fosgrp(at(j)%nfosgrps)%ats(jj)) = at(j)%nfosgrps
      else if (atmres(at(j)%fosgrp(at(j)%nfosgrps)%ats(jj)).eq.(j-1)) then
        which_fosg(at(j)%fosgrp(at(j)%nfosgrps)%ats(jj)) = at(j)%nfosgrps+100
      else
        write(ilog,*) 'Fatal. This is a bug in helper_svgrps_patch(...). Please report.'
        call fexit()
      end if
    end do
  end do
!
  do i=1,nseq
    do j=1,imxfg
      if (allocated(fosbu(i,j)%wts).EQV..true.) then
        deallocate(fosbu(i,j)%wts)
        deallocate(fosbu(i,j)%ats)
      end if
    end do
  end do
  deallocate(fosbu)
  deallocate(grpsprs)
  deallocate(ptchprs)
  deallocate(restlst)
!
end
!
!-------------------------------------------------------------------------------
!

