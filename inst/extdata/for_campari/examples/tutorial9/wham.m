
% Implements the WHAM method, adapted from some initial Andreas-scripts

%Args:
% Main input:  
% A histogram for each window: histo(histo_data1:n,simulation_window)
% histo_ctr: bin centers for histo

function [pmf prob df weights] = wham(histo,histo_ctr,eq_win,ks,T,epsi,refi);

w_bins=length(histo_ctr);
nnodes = length(eq_win);
beta = 4.184/(8.314510*0.001*T);

ubi = zeros(w_bins,nnodes);
prob = zeros(w_bins,1);
pmf = zeros(w_bins,1);
df = zeros(nnodes,1);


hold all;

for i=0:nnodes-1
  %correct for periodicity
  dif=histo_ctr-eq_win(i+1);
  ubi(:,i+1) = ks*(dif).^2.0;
 % plot(histo_ctr,ubi(:,i+1));
 % pause;
end

bl = 1;
bh = w_bins;

%  first construct the initial prob., initializing the Fi to 0
  help1 = zeros(w_bins,1);
  help2 = zeros(nnodes,1);
  for i=1:nnodes
    prob(bl:bh) = prob(bl:bh) + histo(bl:bh,i);
    nlz = sum(histo(bl:bh,i));
    help1(bl:bh) = help1(bl:bh) + nlz*exp(beta*(df(i)-ubi(bl:bh,i)));
  end
  prob(bl:bh) = prob(bl:bh)./help1(bl:bh);
%  second check for consistency in the df
  help2 = df(:);
  for i=1:nnodes
    df(i) = -(1.0/beta)*log ( sum(prob(bl:bh).*exp(-beta*ubi(bl:bh,i))) );
  end
  energy = sqrt( sum((df(:) - help2).^2) );

%  now iterate to convergence
  while (energy > epsi)
    eold = energy;
    help1 = zeros(w_bins,1);
    help2 = zeros(nnodes,1);
    prob(bl:bh) = 0.0;
    for i=1:nnodes
      prob(bl:bh) = prob(bl:bh) + histo(bl:bh,i);
      nlz = sum(histo(bl:bh));
      help1(bl:bh) = help1(bl:bh) + nlz*exp(beta*(df(i)-ubi(bl:bh,i)));
    end
    prob(bl:bh) = prob(bl:bh)./help1(bl:bh);
%   second check for consistency in the df
    help2 = df(:);
    for i=1:nnodes
      df(i) = -(1.0/beta)*log ( sum(prob(bl:bh).*exp(-beta*ubi(bl:bh,i))) );
    end
    energy = sqrt( sum((df(:) - help2).^2) );
  end

% WHAM-estimate
  nzl = df(1);
  df(:) = df(:) - nzl;
  nzl = sum(prob(bl:bh));
  prob(bl:bh) = prob(bl:bh)./nzl;
  pmf(bl:bh) = -(1.0/beta)*log(prob(bl:bh));
  nzl = pmf(refi);
  pmf(bl:bh) = pmf(bl:bh) - nzl;
  
  %Get Normalized weights (make sure we sum to one)

  weights = exp(-beta*(df(:) - df(1)));
  weights = weights ./ sum(weights);

  

