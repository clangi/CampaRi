clear all;

Ts = 298;
zbeq = [0.0 0.1 0.25 0.3 0.4 0.5 0.6 0.75 0.8 0.9 1.0]; %each window in one sim
ks = 150;
nnodes = length(zbeq);
epsi = 0.00005;
bins = 0.005:0.01:0.995;


%Load FB histograms
for i=0:nnodes-1
  fn = sprintf('N_%03d_ZSEC_HIST.dat',i);
  [hdr olap] = hdrload(fn,0);
  hda(:,i+1) = olap(:,3);
end
zbs = bins; %Assume the ZB histo centers NEVER change over any simulation

pmf = zeros(length(zbs));
%Calc the PMF's,Df's.


beta = 4.184/(8.314510*0.001*Ts);
[pmf prob df weights] = wham(hda,zbs,zbeq,ks,Ts,epsi,30);
plot(bins,pmf,'--','linewidth',2);
title('f\beta PMF for Q_{30}');
xlabel('f\beta');
ylabel('PMF(f\beta) in kcal/mol');
axis([0 1 -10 30]);



  

