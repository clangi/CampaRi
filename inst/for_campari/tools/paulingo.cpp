/*

############################################################################
# LICENSE INFO:                                                            #
############################################################################
#    This file is part of CAMPARI.                                         #
#                                                                          #
#    Version 3.0                                                           #
#                                                                          #
#    Copyright (C) 2017, The CAMPARI development team (current and former  #
#                        contributors)                                     #
#                        Andreas Vitalis, Adam Steffen, Rohit Pappu, Hoang #
#                        Tran, Albert Mao, Xiaoling Wang, Jose Pulido,     #
#                        Nicholas Lyle, Nicolas Bloechliger, Marco Bacci,  #
#                        Davide Garolini, Jiri Vymetal                     #
#                                                                          #
#    Website: http://sourceforge.net/projects/campari/                     #
#                                                                          #
#    CAMPARI is free software: you can redistribute it and/or modify       #
#    it under the terms of the GNU General Public License as published by  #
#    the Free Software Foundation, either version 3 of the License, or     #
#    (at your option) any later version.                                   #
#                                                                          #
#    CAMPARI is distributed in the hope that it will be useful,            #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#    GNU General Public License for more details.                          #
#                                                                          #
#    You should have received a copy of the GNU General Public License     #
#    along with CAMPARI.  If not, see <http://www.gnu.org/licenses/>.      #
############################################################################
# AUTHORSHIP INFO:                                                         #
############################################################################
#                                                                          #
# MAIN AUTHOR:   Andreas Vitalis                                           #
#                                                                          #
############################################################################

   a miniature c++-program to compute pauling-derived epsilon parameters
   the poli-vector holds the atomic polarizabilities for the atom types, the radi-vector
   holds the contact distances.
   the program combines the sig_ii to get the sig_ij arithmetically, and then computes eps_ij as:
   eps_ij = constant * pol_i * pol_j / sig_ij^6
   output is directly suitable for use in an ABSINTH parameter file
*/


#include<iostream>
#include<fstream>
#include<cmath>

using namespace std;


int main()
{
   double poli[17]={0.90,1.03,1.34,0.4,0.85,0.69,3.0,0.97,0.4,0.9,0.9,0.18,3.69,0.69,0.4,0.69,0.4};
   double radi[17]={2.70,3.30,3.0,2.0,2.7,3.0,3.6,2.7,2.0,3.2,2.7,1.98,3.62,3.165555296,2.00,3.150656111,2.00};
   double conni=334.60803,val;
   int i,j,naty=17;

   for (i=0;i<naty;i++)
   {
      for (j=i;j<naty;j++)
      {
         val = 0.25*conni*poli[i]*poli[j]/pow(0.5*(radi[i]+radi[j]),6);
         printf("interact         %2d   %2d   %8.6f\n",i+1,j+1,val);
      }
   }
}
