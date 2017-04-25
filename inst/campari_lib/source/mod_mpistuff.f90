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
!                        Nicholas Lyle, Nicolas Bloechliger                !
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
module mpistuff
!
  integer MAXREDIMS
  parameter (MAXREDIMS=35)
!
  type t_commlst
    integer n_s,n_r,n_s2,n_r2
    integer, ALLOCATABLE:: sends(:),recvs(:),sends2(:),recvs2(:)
  end type t_commlst
!
  type(t_commlst) mpi_granullst,mpi_globlst
  logical use_MPIcolls,use_MPIMultiSeed
  integer re_tryswap,re_nbmode,re_freq,re_conditions,re_conddim
  integer re_types(MAXREDIMS),myrank,mpi_nodes,nreolap,re_olcalc,re_aux(10)
  character(MAXSTRLEN) re_infile,re_traceinfile
  integer, ALLOCATABLE:: re_trans(:,:)
  RTYPE, ALLOCATABLE:: re_probs(:,:),re_mat(:,:),re_olap(:,:)
  logical use_REMC,use_MPIAVG,noTI,Tonly,Tthere,reol_all,force_rexyz,force_singlexyz,inst_retr
  integer mpi_cnt_en,mpi_cnt_tor,mpi_cnt_pol,mpi_cnt_xyz,mpi_cnt_sav
  integer mpi_cnt_trcv,mpi_cnt_ens,Tdim,inst_reol,re_velmode,mpi_granularity
  integer, ALLOCATABLE:: mpi_lmap(:,:)
  RTYPE, ALLOCATABLE:: mpi_tvec(:,:),mpi_rbvec(:,:),mpi_cdsvec(:,:,:),mpi_dvec(:,:,:)
  RTYPE, ALLOCATABLE:: mpi_evec(:),mpi_fvec(:),mpi_emat(:,:)
!
end module mpistuff
!

