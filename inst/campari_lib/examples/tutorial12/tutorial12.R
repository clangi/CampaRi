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

# load the data
dat <- read.table("CONTACTMAP.dat",skip=1)

cmap <- array(0.0,dim(dat))
for (i in 1:(dim(dat)[1])) {cmap[i,(1:i)] <- as.double(dat[i,(1:i)]); }

X11(width=5.,height=5.);
image(axes=FALSE,cmap,x=1:(dim(dat)[1]),y=1:(dim(dat)[1]),zlim=c(0,1),col=heat.colors(500)[seq(from=500,to=100,by=-1)],xlab='',ylab='')

box()
axis(side=1,tcl=-0.1,at=1:(dim(dat)[2]),labels=array("",c(dim(dat)[2])))
axis(side=2,tcl=-0.1,at=1:(dim(dat)[2]),labels=array("",c(dim(dat)[2])))

axis(side=1,tcl=0.,at=c(1,4,8,13.5,23.5),labels=c('EAM','HLY','AMP','K+','CL-'),cex.axis=0.6,padj=-2.)
axis(side=2,tcl=0.,at=c(1,4,8,13.5,23.5),labels=c('EAM','HLY','AMP','K+','CL-'),cex.axis=0.6,padj=1.6)

arrows(x0=c(-1,-1,-1,-1),x1=c(100,100,100,100),y0=c(1.5,6.5,8.5,18.5),y1=c(1.5,6.5,8.5,18.5),lty=3)
arrows(y0=c(-1,-1,-1,-1),y1=c(100,100,100,100),x0=c(1.5,6.5,8.5,18.5),x1=c(1.5,6.5,8.5,18.5),lty=3)

legend(x=1,y=(dim(dat)[2]),legend=c("0.00","0.25","0.50","0.75","1.00"),fill=heat.colors(500)[c(500,400,300,200,100)],bg='white')

