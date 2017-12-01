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

reps = 16;
nseq = 21;
minl = 2;
example_directory <- system("pwd",intern=TRUE,ignore.stderr=TRUE,wait=TRUE,input=NULL)

Tsched <- c(250.0,260.0,270.0,280.0,290.0,300.0,310.0,320.0,330.0,340.0,350.0,360.0,370.0,380.0,400.0,420.0)

# load data
if (doread == TRUE) {
  fn <- sprintf("%s/N_000_BB_SEGMENTS_NORM.dat",example_directory);
  dat <- read.table(fn,header=FALSE);
  a <- dim(dat);
  bbcols = a[2];
  bbs <- array(0.0,c(reps,nseq,bbcols));
  fn <- sprintf("%s/N_000_DSSP_NORM.dat",example_directory);
  dat <- read.table(fn,header=FALSE);
  a <- dim(dat);
  dsspcols = a[2];
  dssps <- array(0.0,c(reps,nseq,dsspcols));

  for (rep in 1:reps) {
    fn <- sprintf("%s/N_%03d_BB_SEGMENTS_NORM.dat",example_directory,rep-1);
    dat <- read.table(fn,header=FALSE);
    a <- dim(dat);
    entrs = a[1];
    for (i in 1:entrs) {
      bbs[rep,i,] = as.double(dat[i,]);
    }
  }
  for (rep in 1:reps) {
    fn <- sprintf("%s/N_%03d_DSSP_NORM.dat",example_directory,rep-1);
    dat <- read.table(fn,header=FALSE);
    a <- dim(dat);
    entrs = a[1];
    for (i in 1:entrs) {
      dssps[rep,i,] = as.double(dat[i,]);
    }
  }
  doread <- FALSE;
}

LRs <- array(0.0,c(reps,10))

# compute alpha-content
for (rp1 in 1:reps) {
  for (k in seq(from=nseq,by=-1,to=1)) {
    LRs[rp1,5] = LRs[rp1,5] + k*bbs[rp1,k,4];
    LRs[rp1,6] = LRs[rp1,6] + k*dssps[rp1,k,6];
    if (k >= minl) {
#     the number of H-bonds per segment is not so easily recoverable due to 1) the caps (which may form H-bonds
#     but are not counted), and 2) the definition of what constitutes a helix segment if helical H-bonds are found;
#     here, we assume the maximal number (n-2) rather then the minimal one (n-4), this is also the only consistent
#     setting for DSSP
      if (k > 2) {
        LRs[rp1,3] = LRs[rp1,3] + (k-2.0)*bbs[rp1,k,4]
        LRs[rp1,4] = LRs[rp1,4] + (k-2.0)*dssps[rp1,k,6]
      }
    }
#   this is illegal for DSSP of course for which there are no segments shorter than 4 possible (counts for 2/3 are missed)
    if (k >= minl) {
      LRs[rp1,1] = LRs[rp1,1] + bbs[rp1,k,4]
      LRs[rp1,2] = LRs[rp1,2] + dssps[rp1,k,6]
    }
  }
}

for (rp1 in 1:reps) {
  for (k in 5:6) {LRs[rp1,k] = LRs[rp1,k]/(1.0*nseq);}
}

# compute LR parameters via a MC search
mpower = function(MMi,ppi) {
  ppmax = floor(log2(ppi));
  MMo <- MMi; MMt <- MMi;
  if (ppmax > 0) {
    for (k in 1:ppmax) {MMi <- MMt%*%MMt; MMt <- MMi;}
    if (ppi > (2^ppmax)) {
      for (k in 1:(ppi-(2^ppmax))) {MMi <- MMi%*%MMo;}
    }
  }
  MMo <- MMi
  return(MMo);
}

nb_from_LR = function(v,w,Nr) {
# num. der.
  dw = 1.0e-8
  coremat <- array(0.0,c(3,3))
  coremat2 <- array(0.0,c(3,3))
  dwu = exp(log(w) + dw) - w;
  dwl = exp(log(w) - dw) - w;

  coremat[1,] = c(w+dwu,v,0); coremat[2,] = c(0,0,1); coremat[3,] = c(v,v,1);
  coremat2[1,] = c(w+dwl,v,0); coremat2[2,] = c(0,0,1); coremat2[3,] = c(v,v,1);
  term1 = c(0,0,1) %*% mpower(coremat,Nr) %*% c(0,1,1);
  term2 = c(0,0,1) %*% mpower(coremat2,Nr) %*% c(0,1,1);
  dZm2 = (log(term1) - log(term2))/(2.0*dw);
  return(dZm2)
}

nh_from_LR = function(v,w,Nr) {
# num. der.
  dv = 1.0e-8
  coremat <- array(0.0,c(3,3))
  coremat2 <- array(0.0,c(3,3))
  dvu = exp(log(v) + dv) - v;
  dvl = exp(log(v) - dv) - v;

  coremat[1,] = c(w,v+dvu,0); coremat[2,] = c(0,0,1); coremat[3,] = c(v,v,1);
  coremat2[1,] = c(w,v+dvl,0); coremat2[2,] = c(0,0,1); coremat2[3,] = c(v,v,1);
  term1 = c(0,0,1) %*% mpower(coremat,Nr) %*% c(0,1,1);
  term2 = c(0,0,1) %*% mpower(coremat2,Nr) %*% c(0,1,1);
  dZm2 = (log(term1) - log(term2))/(2.0*dv);
  return(dZm2)
}

v_of_sig = function(sig) {
  coeffs <- c(sig,4*sig,6*sig-1.0,4*sig,sig)
  rts <- polyroot(coeffs)
  rts^2.0/((1.0+rts)^4.0)
  return(rts);
}

rmsd1 <- array(0.0,c(reps,3))
rmsd2 <- array(0.0,c(reps,3))
rmsd1[,3] = nseq;
rmsd2[,3] = nseq;


if (docalc == TRUE) {

  MCsteps <- 2000;
  wssz <- 0.02;
  vssz <- 0.02;
  coremat <- array(0.0,c(3,3))
  coremat2 <- array(0.0,c(3,3))
  vlo = 1.0e-9;
  vhi = 0.5;
  wlo = 1.0e-9;
  whi = 3.5;


  set.seed(system("date +%s",intern=TRUE,ignore.stderr=TRUE,wait=TRUE,input=NULL))

  for (rp1 in 1:reps) {

#   if fresh, use a starting guess, otherwise use the one from the previous temperature
    if (rp1 == 1) {
      wo <- 2.5 #rmsd1[rp1,1];
      vo <- 0.2 #rmsd1[rp1,2];
    }
    else {
      wo <- rmsd1[rp1-1,1];
      vo <- rmsd1[rp1-1,2];
    }
    for (i in 1:MCsteps) {
      if (runif(1,0,1) < 0.1) {wn = runif(1,wlo,whi);}
      else {wn = max(1.0e-9,wo + runif(1,-wssz,wssz));}
      if (runif(1,0,1) < 0.1) {vn = runif(1,vlo,vhi);}
      else {vn = max(1.0e-9,vo + runif(1,-vssz,vssz));}
#     finite difference estimate of partial with w
      dZm1 = nb_from_LR(vn,wn,nseq)
#     finite difference estimate of partial with v12
      dZm2 = nh_from_LR(vn,wn,nseq)
#      writeLines(sprintf("w: %12.5g f: %12.5g <Nh>: %12.5g <Ns>: %12.5g\n",wn,fn,dZm1,dZm2));
#     rmsd w.r.t. to actual data
      bla1 <- sqrt(((dZm1-LRs[rp1,3])/LRs[rp1,3])^2.0 + ((dZm2-LRs[rp1,1])/LRs[rp1,1])^2.0)
      if (bla1 < rmsd1[rp1,3]) {
        LRs[rp1,7] = dZm1
        LRs[rp1,9] = dZm2
        writeLines(sprintf("RB %4d - v: %12.5g w: %12.5g RMSD: %12.5g\n",rp1,vn,wn,bla1))
        rmsd1[rp1,1] = wn;
        wo = wn;
        rmsd1[rp1,2] = vn;
        vo = vn;
        rmsd1[rp1,3] = bla1;
      }
    }
#   the same for DSSP: use BBSEG solution as starting guess
    wo <- rmsd1[rp1,1];
    vo <- rmsd1[rp1,2];
    for (i in 1:MCsteps) {
      if (runif(1,0,1) < 0.1) {wn = runif(1,wlo,whi);}
      else {wn = max(1.0e-9,wo + runif(1,-wssz,wssz));}
      if (runif(1,0,1) < 0.1) {vn = runif(1,vlo,vhi);}
      else {vn = max(1.0e-9,vo + runif(1,-vssz,vssz));}
#     finite difference estimate of partial with w
      dZm1 = nb_from_LR(vn,wn,nseq)
#     finite difference estimate of partial with v12
      dZm2 = nh_from_LR(vn,wn,nseq)
#      writeLines(sprintf("w: %12.5g f: %12.5g <Nh>: %12.5g <Ns>: %12.5g\n",wn,fn,dZm1,dZm2));
#     rmsd w.r.t. to actual data
      bla1 <- sqrt(((dZm1-LRs[rp1,4])/LRs[rp1,4])^2.0 + ((dZm2-LRs[rp1,2])/LRs[rp1,2])^2.0)
      if (bla1 < rmsd2[rp1,3]) {
        LRs[rp1,8] = dZm1
        LRs[rp1,10] = dZm2
        writeLines(sprintf("RB %4d - v: %12.5g w: %12.5g RMSD: %12.5g\n",rp1,vn,wn,bla1))
        rmsd2[rp1,1] = wn;
        wo = wn;
        rmsd2[rp1,2] = vn;
        vo = vn;
        rmsd2[rp1,3] = bla1;
      }
    }
  }
  docalc <- FALSE
}



X11(width=15,height=10);

lnf <- layout(matrix(c(1,2,3,4,5,6),2,3,byrow=TRUE),widths=c(5,5,5),heights=c(5,5))
plot(Tsched,LRs[,5],type='b',col='blue',xlab='Temperature in K',ylab='Net Fractional Helicity',ylim=c(0.0,1.0))
lines(Tsched,LRs[,6],type='b',col='red')
legend(c("BB_SEGMENTS.dat","DSSP.dat","LR Fit to BB_SEGMENTS.dat","LR Fit to DSSP.dat"),x=350.0,y=0.99,lty=c(1,1,2,2),pch=c(1,1,-1,-1),col=c("blue","red","cyan","magenta"))
plot(Tsched,LRs[,3],type='b',col='blue',xlab='Temperature in K',ylab='Mean Number of alpha H-Bonds',ylim=c(0.0,nseq))
lines(Tsched,LRs[,4],type='b',col='red')
lines(Tsched,LRs[,7],type='l',col='cyan',lty=2)
lines(Tsched,LRs[,8],type='l',col='magenta',lty=2)

plot(Tsched,LRs[,1],type='b',col='blue',xlab='Temperature in K',ylab='Mean Number of a-Helical Segments',ylim=c(0.0,max(LRs[,1],LRs[,2])))
lines(Tsched,LRs[,2],type='b',col='red')
lines(Tsched,LRs[,9],type='l',col='cyan',lty=2)
lines(Tsched,LRs[,10],type='l',col='magenta',lty=2)

plot(Tsched,rmsd1[,1],type='b',col='blue',xlab='Temperature in K',ylab='Lifson-Roig Propagation Parameter',ylim=c(0.0,max(rmsd1[,1],rmsd2[,1])))
lines(Tsched,rmsd2[,1],type='b',col='red')
plot(Tsched,rmsd1[,2],type='b',col='blue',xlab='Temperature in K',ylab='Lifson-Roig Nucleation Parameter',ylim=c(0.0,max(rmsd1[,2],rmsd2[,2])))
lines(Tsched,rmsd2[,2],type='b',col='red')




