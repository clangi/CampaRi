%--------------------------------------------------------------------------%
% LICENSE INFO:                                                            %
%--------------------------------------------------------------------------%
%    This file is part of CAMPARI.                                         %
%                                                                          %
%    Version 2.0                                                           %
%                                                                          %
%    Copyright (C) 2014, The CAMPARI development team (current and former  %
%                        contributors)                                     %
%                        Andreas Vitalis, Adam Steffen, Rohit Pappu, Hoang %
%                        Tran, Albert Mao, Xiaoling Wang, Jose Pulido,     %
%                        Nicholas Lyle, Nicolas Bloechliger                %
%                                                                          %
%    Website: http://sourceforge.net/projects/campari/                     %
%                                                                          %
%    CAMPARI is free software: you can redistribute it and/or modify       %
%    it under the terms of the GNU General Public License as published by  %
%    the Free Software Foundation, either version 3 of the License, or     %
%    (at your option) any later version.                                   %
%                                                                          %
%    CAMPARI is distributed in the hope that it will be useful,            %
%    but WITHOUT ANY WARRANTY; without even the implied warranty of        %
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         %
%    GNU General Public License for more details.                          %
%                                                                          %
%    You should have received a copy of the GNU General Public License     %
%    along with CAMPARI.  If not, see <http://www.gnu.org/licenses/>.      %
%--------------------------------------------------------------------------%
% AUTHORSHIP INFO:                                                         %
%--------------------------------------------------------------------------%
%                                                                          %
% MAIN AUTHOR:   Andreas Vitalis                                           %
%                                                                          %
%--------------------------------------------------------------------------%

function [conc, dens] = box2conc(Np,Mass,rv)
%
%BOX2CONC: A Simple Conversion Tool to Get Concentration/Density from Box Parameters
%
%   Usage: box2conc(N_particles,Mass_tot,box_vec)
%
%   N_particles is the number of particles (scalar), Mass is the total mass (scalar)
%   in g/mol and box_vec is either a scalar (sphere radius), or a 3D-vector (rectangular
%   box side lengths) given in Angstroms.
%
%   Negative values are not allowed for particle number or mass, neither are resultant
%   volumes less or equal to zero (sign of box vector is not checked).
%   
conc = -1.0;
dens = -1.0;
[a b] = size(rv);

if (a ~= 1) && (b ~= 1) return; end
if (a == 3) && (b == 1)
  cub = 1;
elseif (a == 1) && (b == 3)
  cub = 1; rv = rv';
elseif (a == 1) && (b == 1)
  cub = 0;
else
  box2conc=[-1 -1]; return;
end

[a b] = size(Mass);
if (a ~= 1) || (b ~= 1) return; end
if (Mass < 0.0) return; end

[a b] = size(Np);
if (a ~= 1) || (b ~= 1) return; end
if (Np < 0.0) return; end


clear Na; clear Atol; clear voli;
Na = 6.0221367e23;
Atol = 1.0e-27;

if (cub == 1)
  voli = Atol*rv(1)*rv(2)*rv(3);
else
  voli = Atol*(4./3.)*pi*(rv^3);
end

if (voli <= 0.0); return; end

conc = Np/(Na*voli);
dens = 0.001*Mass/(Na*voli);

return;


