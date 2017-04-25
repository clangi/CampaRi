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
! acc             : acceptance counts type
! naccept         : fyc move acceptance counts
! nchi            : chi move acceptance counts
! nfy             : pivot move acceptance counts
! nre             : replica exchange acceptance counts
! nomega          : omega move acceptance counts
! nlct            : linear combination of torsions move acceptance counts
! nrb             : rigid body move acceptance counts
! nnuc            : nucleotide move acceptance counts
! ninsert         : particle insertion move acceptance counts
! ndelete         : particle deletion move acceptance counts
! nidentitychange : identity swap move acceptance counts
! npucker         : sugar pucker move acceptance counts
! fy(:)           : backbone move acceptance counts per residue
! chi(:)          : chi move acceptance counts per residue
! nuc(:)          : nucleotide move accpetance counts per residue
! pucker(:)       : pucker move acceptance counts per residue
! rigid(:)        : rigid move acceptance counts per residue
! omega(:)        : omega move acceptance counts per residue
! cr(:)           : concerted rotation move acceptance counts per residue
! insert(:)       : insertion move acceptance counts per molecule
! delete(:)       : deletion move acceptance counts per molecule
! permute(:)      : identitfy swap move acceptance counts per molecule
!
module accept
!
  type t_accept
    integer naccept,nchi,nfy,nre,nomega,nlct,nrb,nnuc,npucker,nnuccr
    integer ninsert,ndelete,nidentitychange,nujcr,ndjcr,ndocr,nsjcr,nfyc
    integer nclurb,ntrans,nrot,nph,nother
    integer, ALLOCATABLE:: fy(:),chi(:),nuc(:)
    integer, ALLOCATABLE:: rigid(:),omega(:),cr(:),other(:)
    integer, ALLOCATABLE:: insert(:),delete(:),permute(:),pucker(:)
  end type t_accept
  type(t_accept):: acc
!
end module accept
!
