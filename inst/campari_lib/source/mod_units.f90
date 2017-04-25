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
! MAIN AUTHOR:   Rohit Pappu                                               !
! CONTRIBUTIONS: Andreas Vitalis, Albert Mao                               !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
!
! avogadro      : Na in # * mol^-1
! gasconst      : R in kcal * mol^-1 * K^-1
! electric      : standard Coulomb interaction between unit charges at 1A separation in kcal * mol^-1
! debye         : standard dipole moment of two unit charges at 1A separation in Debye
! convdens      : conversion from a.m.u * Ang^-3 to g * cm^-3
! u_dyn_kb      : kb in  g * A^2 * ps^-2 * K^-1 * mol^-1 
! u_dyn_fconv   : conversion from kcal * mol^-1 * A-1 to g * A^1 * ps^-2 * mol^-1
! u_dyn_fconv_r : conversion from kcal * mol^-1 * rad-1 to g * A^2 * deg^-1 * ps^-2 * mol^-1
! u_dyn_virconv : conversion from kcal * mol^-1 * A-1 to bar * A^2
! planckconstant: h in Angstrom sqrt(gram * kcal) / mol

!
module units
!
  RTYPE avogadro,gasconst,electric,debye,convdens,planckconstant
  RTYPE u_dyn_kb,u_dyn_fconv,u_dyn_fconv_r,u_dyn_virconv
! 
  parameter (avogadro=6.0221418d+23)
  parameter (gasconst=1.98717623d-3)
  parameter (electric=332.05382d0)
  parameter (debye=4.8033324d0)
  parameter (convdens=1.66053886d0)
  parameter (u_dyn_kb=0.8314510d0)
  parameter (u_dyn_fconv=418.4d0)
  parameter (u_dyn_fconv_r=23972.554d0)
  parameter (u_dyn_virconv=6.9477001d+4)
  parameter (planckconstant=1.9501384275311227d0)
!
end module units
!

