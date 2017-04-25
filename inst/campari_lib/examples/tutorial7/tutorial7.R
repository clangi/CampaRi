############################################################################
# LICENSE INFO:                                                            #
############################################################################
#    This file is part of CAMPARI.                                         #
#                                                                          #
#    Version 2.0                                                           #
#                                                                          #
#    Copyright (C) 2014, The CAMPARI development team (current and former  #
#                        contributors)                                     #
#                        Andreas Vitalis, Adam Steffen, Rohit Pappu, Hoang #
#                        Tran, Albert Mao, Xiaoling Wang, Jose Pulido,     #
#                        Nicholas Lyle, Nicolas Bloechliger                #
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

reps = 24;
kbT = 0.298*8.314510/4.184;
beta = 1.0/kbT;
df <- array(0.0,c(reps,reps));
dfsc <- array(0.0,c(reps,reps));

ljlr <- -0.0619;

example_directory <- system("pwd",intern=TRUE,ignore.stderr=TRUE,wait=TRUE,input=NULL)

sched <- c(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.65,0.7,0.75,0.8,0.85,0.875,0.9,0.925,0.95,0.975,1.0,1.25,1.5,1.65,1.8,1.9,2.0)
schedm0 <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.65,0.7,0.75,0.8,0.85,0.875,0.9,0.925,0.95,0.975,1.0,1.25,1.5,1.65,1.8,1.9,2.0)


# load data
if (doread == TRUE) {
  fn <- sprintf("%s/N_001_EVEC.dat",example_directory);
  dat <- read.table(fn);
  a <- dim(dat);
  entrs = a[1];

  ev <- array(0.0,c(reps,entrs,reps));
  for (rep in 1:reps) {
    fn <- sprintf("%s/N_%03d_EVEC.dat",example_directory,rep-1);
    dat <- read.table(fn);
    for (i in 1:entrs) {
      if ((rep > 1) && (rep < reps)) {for (k in (rep-1):(rep+1)) {ev[rep,i,k] = as.double(dat[i,(k-(rep-1)+1)]);}}
      if (rep == reps) {ev[rep,i,rep] = as.double(dat[i,2]); ev[rep,i,rep-1] = as.double(dat[i,1]);}
      if (rep == 1) {ev[rep,i,rep] = as.double(dat[i,1]); ev[rep,i,rep+1] = as.double(dat[i,2]);}
    }
  }

  fn <- sprintf("%s/SELFCORR/N_001_EVEC.dat",example_directory);
  dat <- read.table(fn);
  a <- dim(dat);
  entrs = a[1];

  evsc <- array(0.0,c(reps,entrs,reps));
  for (rep in 1:reps) {
    fn <- sprintf("%s/SELFCORR/N_%03d_EVEC.dat",example_directory,rep-1);
    dat <- read.table(fn);
    for (i in 1:entrs) {
      if ((rep > 1) && (rep < reps)) {for (k in (rep-1):(rep+1)) {evsc[rep,i,k] = as.double(dat[i,(k-(rep-1)+1)]);}}
      if (rep == reps) {evsc[rep,i,rep] = as.double(dat[i,2]); evsc[rep,i,rep-1] = as.double(dat[i,1]);}
      if (rep == 1) {evsc[rep,i,rep] = as.double(dat[i,1]); evsc[rep,i,rep+1] = as.double(dat[i,2]);}
    }
  }
  doread <- FALSE;
}

# BAR
for (rp1 in 1:(reps-1)) {
  rp2 = rp1+1;
  eps = 0.2;
  baro = 0.0;
  while (abs(eps) > 0.000005) {
    barn = baro + kbT*log( sum(1.0/(1.0 + exp(beta*(ev[rp2,,rp1]-ev[rp2,,rp2]+baro)))) / sum((1.0/(1.0 + exp(beta*(ev[rp1,,rp2]-ev[rp1,,rp1]-baro))))) );
    eps = barn - baro;
    baro = barn;
  }
  df[rp1,rp2] = baro;

  eps = 0.2;
  baro = 0.0;
  while (abs(eps) > 0.000005) {
    barn = baro + kbT*log( sum(1.0/(1.0 + exp(beta*(ev[rp1,,rp2]-ev[rp1,,rp1]+baro)))) / sum((1.0/(1.0 + exp(beta*(ev[rp2,,rp1]-ev[rp2,,rp2]-baro))))) );
    eps = barn - baro;
    baro = barn;
  }
  df[rp2,rp1] = baro;
}

for (rp1 in 1:(reps-1)) {
  rp2 = rp1+1;
  eps = 0.2;
  baro = 0.0;
  while (abs(eps) > 0.000005) {
    barn = baro + kbT*log( sum(1.0/(1.0 + exp(beta*(evsc[rp2,,rp1]-evsc[rp2,,rp2]+baro)))) / sum((1.0/(1.0 + exp(beta*(evsc[rp1,,rp2]-evsc[rp1,,rp1]-baro))))) );
    eps = barn - baro;
    baro = barn;
  }
  dfsc[rp1,rp2] = baro;

  eps = 0.2;
  baro = 0.0;
  while (abs(eps) > 0.000005) {
    barn = baro + kbT*log( sum(1.0/(1.0 + exp(beta*(evsc[rp1,,rp2]-evsc[rp1,,rp1]+baro)))) / sum((1.0/(1.0 + exp(beta*(evsc[rp2,,rp1]-evsc[rp2,,rp2]-baro))))) );
    eps = barn - baro;
    baro = barn;
  }
  dfsc[rp2,rp1] = baro;
}


dfraw <- 0.0
selfcorrs <- 0.0
for (rp1 in 2:reps) {
  dfraw = dfraw + df[rp1-1,rp1]
  selfcorrs = selfcorrs + dfsc[rp1-1,rp1]
}
dftot = dfraw + ljlr - selfcorrs

writeLines(sprintf("The final BAR estimate including all corrections is %g kcal/mol.",dftot));
writeLines(sprintf("This includes long-range LJ corrections contributing %g kcal/mol and self-corrections of %g kcal/mol (being subtracted out).",ljlr,selfcorrs));
writeLines(sprintf("The uncorrected FES estimate via BAR is %g kcal/mol.",dfraw));




