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
! CONTRIBUTIONS: Nicolas Bloechliger                                       !
! WRAPPER: Davide Garolini                                                 !
!                                                                          !
!--------------------------------------------------------------------------!

! -----------------------------------------------------------------------------
! This function aims to reduce the fringe regions by collapsing the tree leaves
! on their branches. In this way we have less grouping of fringe regions
!
subroutine contract_mst_fortran(n_snaps,mnb,alnbs,alst,aldis,nrnds,istats)
!
  implicit none
!
  integer, INTENT(IN) :: n_snaps,mnb,nrnds,alnbs(n_snaps),alst(n_snaps,mnb)
  real, INTENT(INOUT) :: aldis(n_snaps,mnb)
  integer, INTENT(OUT) :: istats(2)
!
  integer i, j, kk, thej, ll
!
  logical, ALLOCATABLE :: terminal(:)
!
  allocate(terminal(n_snaps))
  terminal(:) = .false.
!
  istats(2) = 0
  do i=1,n_snaps
    if (alnbs(i).eq.1) then
      terminal(i) = .true.
      istats(2) = istats(2) + 1
    end if
  end do
  ! write(*,*) nrnds,n_snaps,istats(2)
!
  istats(1) = 0
  do kk=1,nrnds
    do i=1,n_snaps
      if (terminal(i).EQV..true.) then
        do j=1,alnbs(i)
          if (aldis(i,j).ge.0.0) then
            thej = alst(i,j)
            if (abs(aldis(i,j)).le. 10.0E-12) then
              aldis(i,j) = -10.0*TINY(aldis(i,j))
            else
              aldis(i,j) = -1.0/aldis(i,j)
            end if
            istats(1) = istats(1) + 1
            exit
          end if
        end do
        do ll=1,alnbs(thej)
          if (alst(thej,ll).eq.i) then
            aldis(thej,ll) = aldis(i,j)
            exit
          end if
        end do
      end if
    end do
!   reset terminal status
    terminal(:) = .false.
    istats(2) = 0
    do i=1,n_snaps
      thej = 0
      do j=1,alnbs(i)
        if (aldis(i,j).ge.0.0) thej = thej + 1
      end do
      if (thej.eq.1) then
        terminal(i) = .true.
        istats(2) = istats(2) + 1
      end if
    end do
    if (istats(2).eq.0) exit
  end do
!
  deallocate(terminal)
!
end
