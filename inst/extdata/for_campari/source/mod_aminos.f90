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
! CONTRIBUTIONS: Rohit Pappu                                               !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
! integer codes for possible N-terminal residues: 1-22,24,27-28,33-35,51,64-87,108-119
! integer codes for possible C-terminal residues: 1-22,24,29-30,33-35,51,64-75,108-119
! integer codes for mandatory single residue mol.s: 36-50,52-63,88-107
!
!
module aminos
!
  integer MAXAMINO
  parameter (MAXAMINO=119)
!
  character(1) aminopolty(MAXAMINO)
  character(3) amino(MAXAMINO)
  data amino  / 'GLY','ALA','VAL','LEU','ILE','SER','THR','CYS',&
               &'PRO','PHE','TYR','TRP','HIP','HID','HIE','ASP',&
               &'ASN','GLU','GLN','MET','LYS','ARG','ORN','AIB',&
               &'PCA','UNK','ACE','FOR','NME','NH2','GAM','HYP',&
               &'ABA','NVA','NLE','NA+','CL-','URE','SPC','T3P',&
               &'NMF','NMA','ACA','PPA','T4P','CH4','FOA','DMA',&
               &'MOH','PCR','DAB','K+ ','BR-','CS+','I- ','NH4',&
               &'AC-','GDN','CYT','URA','THY','ADE','GUA','R5P',&
               &'D5P','RPU','RPC','RPT','RPA','RPG','DPU','DPC',&
               &'DPT','DPA','DPG','RIB','DIB','RIU','RIC','RIT',&
               &'RIA','RIG','DIU','DIC','DIT','DIA','DIG','PRP',&
               &'NBU','IBU','TOL','EOH','MSH','EMT','IMD','IME',&
               &'MIN','1MN','2MN','LCP','NO3','T5P','T4E','BEN',&
               &'NAP','PUR','O2 ','GLH','ASH','CYX','TYO','LYD',&
               &'PTR','SEP','TPO','KAC','KM1','KM2','KM3'/
  data aminopolty / 'P','P','P','P','P','P','P','P','P','P','P','P',&
                 &'P','P','P','P','P','P','P','P','P','P','P','P',&
                 &'P','X','X','X','X','X','P','P','P','P','P','X',&
                 &'X','X','X','X','X','X','X','X','X','X','X','X',&
                 &'X','X','P','X','X','X','X','X','X','X','X','X',&
                 &'X','X','X','N','N','N','N','N','N','N','N','N',&
                 &'N','N','N','N','N','N','N','N','N','N','N','N',&
                 &'N','N','N','X','X','X','X','X','X','X','X','X',&
                 &'X','X','X','X','X','X','X','X','X','X','X','P',&
                 &'P','P','P','P','P','P','P','P','P','P','P'/
  save amino,aminopolty
!
  RTYPE rescontlen(MAXAMINO)
  data rescontlen &
             &/  3.64, 3.64, 3.64, 3.64, 3.64, 3.64, 3.64, 3.64,&
                &3.64, 3.64, 3.64, 3.64, 3.64, 3.64, 3.64, 3.64,&
                &3.64, 3.64, 3.64, 3.64, 3.64, 3.64, 3.64, 3.64,&
                &3.64, 0.00, 2.80, 1.33, 2.50, 1.33, 3.64, 3.64,&
                &3.64, 3.64, 3.64, 0.00, 0.00, 2.10, 1.10, 1.10,&
                &2.30, 3.40, 2.20, 3.50, 1.10, 1.10, 2.10, 3.50,&
                &2.00, 4.70, 3.64, 0.00, 0.00, 0.00, 0.00, 1.00,&
                &2.30, 2.10, 3.80, 3.80, 4.10, 3.80, 4.30, 7.00,&
                &7.00, 7.00, 7.00, 7.00, 7.00, 7.00, 7.00, 7.00,&
                &7.00, 7.00, 7.00, 4.00, 4.00, 4.00, 4.00, 4.00,&
                &4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00,&
                &4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00,&
                &4.00, 4.00, 4.00, 4.00, 4.00, 1.10, 1.10, 4.00,&
                &4.00, 4.00, 1.10, 3.64, 3.64, 3.64, 3.64, 3.64,&
                &3.64, 3.64, 3.64, 3.64, 3.64, 3.64, 3.64/
  save rescontlen
!
end module aminos
!
