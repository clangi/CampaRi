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
% CONTRIBUTIONS: Adam Steffen                                              %
%                                                                          %
%--------------------------------------------------------------------------%

function exitval(mini,maxi,centro,steepo)
%
%SIGMAINTERPOL: Sigmoidal interpolation between limits
%
%   Usage: sigmainterpol(min,max,center,tau)
%
%   Calculates and plots a sigmoidal function between the limits [min max]
%   (assumed within [0 1]) such as to map it to the interval [0 1] to form
%   a continous (but not smooth) interpolation with variable parameters.
%
%   The center ([0 1]) determines the location of the mid-point of the
%   sigmoidal function relative to the theoretical limits (0.0 corresponds
%   to min, 1.0 corresponds to max).
%
%   Tau determines the inherent decay constant of the exponential:
%   f(x) ~ 1.0 / (1.0 + exp(-x/tau))
%
%   All input parameters are assumed to be scalars.
%   
[a b] = size(mini);
if (a ~= 1) || (b ~= 1) exitval=-1; return; end

[a b] = size(maxi);
if (a ~= 1) || (b ~= 1) exitval=-1; return; end

[a b] = size(centro);
if (a ~= 1) || (b ~= 1) exitval=-1; return; end

[a b] = size(steepo);
if (a ~= 1) || (b ~= 1) exitval=-1; return; end


clear stretchi; clear centri; clear sosi; clear offi;
sosi = zeros(1,100);

centri = centro*maxi + (1.0-centro)*mini;
stretchi = 1.0/ (1.0/(1.0 + exp(-(maxi-centri)/steepo)) - 1.0/(1.0 + exp(-(mini-centri)/steepo)));
offi = 1.0 - (1.0/(1.0 + exp(-(maxi-centri)/steepo)) - 0.5)*stretchi;

for i=1:100
  savi = (i-0.5)*0.01;
  if (savi <= mini) sosi(i) = 0.0;
  elseif (savi >= maxi) sosi(i) = 1.0;
  else
   sosi(i) = 1.0/(1.0 + exp(-(savi-centri)/steepo));
   sosi(i) = (sosi(i)-0.5)*stretchi + offi;
  end
end

plot(0.005:0.01:0.995,sosi,'ro-');

exitval = 0;
return;


