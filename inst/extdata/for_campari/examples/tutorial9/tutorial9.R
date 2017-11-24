
kss <- 150.0 # force constant (FMCSC_ZS_FR_KB)
scalef <- 1.0 # outside scale factor (FMCSC_SC_ZSEC)

zbeq = seq(from=0.0,to=1.0,by=0.1); # choice of window centers
zbs = seq(from=0.005,to=0.995,by=0.01); # default binning in N_*_ZSEC_HIST.dat

beta = 4.184/(8.314510*0.298); # calculate inverse temperature from FMCSC_TEMP

epsi = 0.00001; # convergence criterion
bins = length(zbs)
numnodes = length(zbeq)


zsecs <- array(0,c(numnodes,bins))
for (k in 1:numnodes) {
  tmps <- paste(sprintf("N_%03d_ZSEC_HIST.dat",k-1))
  d=read.table(tmps,skip=1);
  zsecs[k,] = as.double(d[,3])
}


ubi <- array(0.0,c(numnodes,bins))
prob <- array(0.0,c(bins))
pmf <- array(0.0,c(bins))
df <- array(0.0,c(numnodes))

for (k in 1:numnodes) {
  ubi[k,] = kss*((zbs[]-zbeq[k])^2.0); # the bias potential
}


bl = 1
bh = 100
hadsome <- FALSE
for (k in 1:bins) {
  if (sum(zsecs[,k]) > 0.0) {hadsome = TRUE}
  if ((sum(zsecs[,k]) <= 0.0) && (hadsome == FALSE)) {
    bl = k + 1
  }
  if ((sum(zsecs[,k]) <= 0.0) && (hadsome == TRUE)) {
    bh = k - 1
    break
  }
}
if (bl >= bh) {print("Fatal"); return}
help1 <- array(0.0,bins)
help2 <- array(0.0,max(numnodes))
for (k in 1:numnodes) {
  prob[bl:bh] = prob[bl:bh] + zsecs[k,bl:bh]
  nlz <- sum(zsecs[k,bl:bh])
  help1[bl:bh] = help1[bl:bh] + nlz*exp(beta*(df[k] - ubi[k,bl:bh]))
}
prob[bl:bh] = prob[bl:bh]/help1[bl:bh]
## check consistency 
help2 = df[]
for (k in 1:numnodes) {
  df[k] = -(1.0/beta)*log ( sum(prob[bl:bh] * exp(-beta*ubi[k,bl:bh])) )
}
energy <- sqrt( sum( (df[] - help2)^2.0) )
## iterate
niters <- 0
while (energy > epsi) {
  niters <- niters + 1
  eold = energy
  help1 <- array(0.0,bins)
  help2 <- array(0.0,max(numnodes))
  prob[bl:bh] = 0.0
  for (k in 1:numnodes) {
    prob[bl:bh] = prob[bl:bh] + zsecs[k,bl:bh]
    nlz <- sum(zsecs[k,bl:bh])
    help1[bl:bh] = help1[bl:bh] + nlz*exp(beta*(df[k] - ubi[k,bl:bh]))
  }
  prob[bl:bh] = prob[bl:bh]/help1[bl:bh]
  help2 = df[]
  for (k in 1:numnodes) {
    df[k] = -(1.0/beta)*log ( sum(prob[bl:bh] * exp(-beta*ubi[k,bl:bh])) )
  }
  energy <- sqrt( sum( (df[] - help2)^2.0) )
}
nzl = df[1]
df[] = df[] - nzl;
nzl = sum(prob[bl:bh]);
prob[bl:bh] = prob[bl:bh]/nzl
pmf[bl:bh] = -(1.0/beta)*log( prob[bl:bh] )
nzl = min(pmf[bl:bh])
pmf[] = pmf[] - nzl

dev.new()
plot(x=zbs,y=zsecs[1,],xlab=expression(paste(sep="",f[beta])),ylab="Frequency",col=heat.colors(2*numnodes)[1],ylim=c(0.0,max(zsecs)),type='l',lwd=2.,main="Biased distributions")
for (k in 2:numnodes) {
  lines(x=zbs,y=zsecs[k,],col=heat.colors(2*numnodes)[k],lwd=2.);
}

dev.new()
plot(x=zbs[bl:bh],y=pmf[bl:bh],xlim=c(0,1),xlab=expression(paste(sep="",f[beta])),ylab="PMF in kcal/mol",col='black',ylim=c(min(pmf),max(pmf)),type='l',lwd=2.,main="Potential of mean force")
lines(type='s',y=c(df,df[numnodes])-min(df),x=c(zbeq-0.5*(zbeq[2]-zbeq[1]),zbeq[numnodes]),lty=2,lwd=2.,col='cyan4')
legend(x='top',legend=c("Potential of mean force","Free energy of windows"),col=c("black","cyan4"),lty=c(1,2),lwd=2,pch=-1)
