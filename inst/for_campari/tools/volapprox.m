%--------------------------------------------------------------------------%
% LICENSE INFO:                                                            %
%--------------------------------------------------------------------------%
%    This file is part of CAMPARI.                                         %
%                                                                          %
%    Version 3.0                                                           %
%                                                                          %
%    Copyright (C) 2017, The CAMPARI development team (current and former  %
%                        contributors)                                     %
%                        Andreas Vitalis, Adam Steffen, Rohit Pappu, Hoang %
%                        Tran, Albert Mao, Xiaoling Wang, Jose Pulido,     %
%                        Nicholas Lyle, Nicolas Bloechliger, Marco Bacci,  %
%                        Davide Garolini, Jiri Vymetal                     %
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
% CONTRIBUTIONS: Adam Steffen                                              %
%                                                                          %
%--------------------------------------------------------------------------%

function slup = exitval(r1,r2,rw,dres)
%
%VOLAPPROX: Comparison of Linear Approximation for Volume Overlap to
%               Analytical Solution
%
%   Usage: volapprox(ra,rb,rw,res)
%
%   Calculates the overlap volume of a sphere or radius rb with the
%   spherical shell of radius rw around a sphere of radius ra as a
%   function of the two spheres. The distance resolution is set by
%   res. The maximum value of relevance is ra+rb+rw.
%   The output matrix has 5 columns, column 1 has the distance values,
%   column 2 has the analytical values, and column 3 the numerical 
%   values. The inverted problem (spherical shell around sphere of
%   radius rb with second sphere of radius ra) is given in columns
%   4 and 5, respectively.
%   

[a b] = size(r1);
if (a ~= 1) || (b ~= 1) exitval=-1; return; end

[a b] = size(r2);
if (a ~= 1) || (b ~= 1) exitval=-1; return; end

[a b] = size(rw);
if (a ~= 1) || (b ~= 1) exitval=-1; return; end


dmax = r1+r2+rw;
d = dres:dres:dmax;

[a nbins] = size(d);

Vol = zeros(nbins,5);
vol1 = (4./3.)*pi*r1^3.0;
vol2 = (4./3.)*pi*r2^3.0;

for i=1:nbins
  if (d(i) > dmax - 2*r2)
    Vol(i,3) = vol2*(dmax - d(i))/(2.0*r2);
    Vol(i,2) = (1.0/(12.0*d(i))) * pi * ((r1+rw+r2-d(i))^2.0) * (d(i)^2.0 + 2.0*d(i)*r2 - 3.0*r2^2.0 + 2*d(i)*(r1+rw) + 6.0*(r1+rw)*r2 - 3.0*(r1+rw)^2.0);
  elseif (d(i) < r1 + r2)
    Vol(i,3) = vol2*d(i)/(r1+r2);
    if (r1 > r2)
      if (d(i) < r1-r2)
        Vol(i,2) = 0.0;
      else
        Vol(i,2) = vol2 - (1.0/(12.0*d(i))) * pi * ((r1+r2-d(i))^2.0) * (d(i)^2.0 + 2.0*d(i)*r2 - 3.0*r2^2.0 + 2*d(i)*r1 + 6.0*r1*r2 - 3.0*r1^2.0);
      end
    else
      if (d(i) < r2-r1)
        Vol(i,2) = vol2 - vol1;
      else
        Vol(i,2) = vol2 - (1.0/(12.0*d(i))) * pi * ((r2+r1-d(i))^2.0) * (d(i)^2.0 + 2.0*d(i)*r1 - 3.0*r1^2.0 + 2*d(i)*r2 + 6.0*r1*r2 - 3.0*r2^2.0);
      end
    end
  else
    Vol(i,3) = vol2;
    Vol(i,2) = vol2;
  end
end

for i=1:nbins
  if (d(i) > dmax - 2*r1)
    Vol(i,5) = vol1*(dmax - d(i))/(2.0*r1);
    Vol(i,4) = (1.0/(12.0*d(i))) * pi * ((r1+rw+r2-d(i))^2.0) * (d(i)^2.0 + 2.0*d(i)*r1 - 3.0*r1^2.0 + 2*d(i)*(r2+rw) + 6.0*(r2+rw)*r1 - 3.0*(r2+rw)^2.0);
  elseif (d(i) < r1 + r2)
    Vol(i,5) = vol1*d(i)/(r1+r2);
    if (r2 > r1)
      if (d(i) < r2-r1)
        Vol(i,4) = 0.0;
      else
        Vol(i,4) = vol1 - (1.0/(12.0*d(i))) * pi * ((r1+r2-d(i))^2.0) * (d(i)^2.0 + 2.0*d(i)*r1 - 3.0*r1^2.0 + 2*d(i)*r2 + 6.0*r1*r2 - 3.0*r2^2.0);
      end
    else
      if (d(i) < r1-r2)
        Vol(i,4) = vol1 - vol2;
      else
        Vol(i,4) = vol1 - (1.0/(12.0*d(i))) * pi * ((r2+r1-d(i))^2.0) * (d(i)^2.0 + 2.0*d(i)*r2 - 3.0*r2^2.0 + 2*d(i)*r1 + 6.0*r1*r2 - 3.0*r1^2.0);
      end
    end
  else
    Vol(i,5) = vol1;
    Vol(i,4) = vol1;
  end
end

Vol(:,1) = d;
slup = Vol;
return;


