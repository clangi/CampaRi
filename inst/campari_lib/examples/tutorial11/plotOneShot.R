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
# MAIN AUTHOR:   Nicolas Bloechliger                                       #
#                                                                          #
############################################################################

doread <- TRUE
if (doread == TRUE) {
  progind <- read.table("PROGIDX_000000000001.dat")[, c(3, 4)]
  fyc <- read.table("NBU.fyc")[, c(2:4)]
  # Determine, based on fyc, which snapshots are eclipsed:
  e <- (abs(fyc) < 30 | abs(abs(fyc) - 120) < 30)
  eclipsed <- e[, 1] | e[, 2] | e[, 3]
  e <- NULL

  doread <- FALSE
}

tmps <- paste(sep="", "oneShot.eps")
setEPS()
postscript(tmps, width=3.54330709, height=2.5, paper="special")
par(mar=c(1.3, 1.4, 0.1, 0.1))
par(mgp=c(0, -0.4, 0))
par(ps=6)

cstored <- dim(progind)[1]  # Number of snapshots in trajectory
xr = range(c(1, cstored))  # Range of x-values in plot
ymin = min(-log(progind[, 2] / cstored))
ymax = max(-log(progind[100:(cstored-100), 2] / cstored))
yr = range(c(ymin, ymax))  # Range of y-values in plot

# Print percentage of eclipsed snapshots:
print(paste(sep="", "Note: ", format(100 * sum(eclipsed) / cstored, digits=2), "% of all snapshots are eclipsed"))  

# Plot one-shot profile:
plot(1:cstored, -log(progind[, 2] / cstored), xlim=xr, ylim=yr, xlab="", ylab="", axes=FALSE, pch=16, cex=0.1, col=(1 + eclipsed[progind[, 1]]))
box()

# Add axis labels:
axis(1, at=c(0, 25000, 50000, 75000, 100000), labels=c(0, expression(2.5%.%10^4), expression(5%.%10^4), expression(7.5%.%10^4), expression(10^5)), tck=.01, padj=0.3)
mtext("Progress Index", side=1, line=0.3, cex=4/3)
yy <- format(c(ymin, ymin + (ymax - ymin) * 0.25, ymin + (ymax - ymin) * 0.5), digits=2)
axis(2, at=yy, las=3, tck=.01, padj=-1)
mtext(expression("ln(("*italic(tau["SA"]+tau["AS"])*")/2)"), side=2, line=0.4, cex=4/3, at=yy[2])

# Add dihedral annotation:
ro.cols <- colorRampPalette(c(rgb(198, 219, 239, maxColorValue=255), rgb(107, 174, 214, maxColorValue=255), rgb(8, 81, 156, maxColorValue=255)), space = "Lab")  # Defines set of colors
xx <- seq(from=1, by=10, to=cstored)  # To reduce image size, only plot dihedral annotation for every 10th snapshot along the progress index
yy <- ymin + (ymax - ymin) * c(2/3, 11/15, 12/15, 13/15)  # Vertical position of dihedral annotation
image(xx, y=yy, z=((fyc[progind[xx, 1], ] %% 360 < 120) + (fyc[progind[xx, 1], ] %% 360 < 240) + 1), col=ro.cols(3), bg="white", axes=FALSE, xlim=xr, add=TRUE)  # Plot dihedral annotation
legend(xpd=TRUE, x=(cstored / 2), y=(ymin + (ymax - ymin) * 14 / 15), xjust=0.5, yjust=0.5, legend=c(expression("gauche"^"-"), "anti", expression( "gauche"^"+")), fill=ro.cols(3), ncol=3, bty="n", cex=4/3, x.intersp=0.4)  # Legend for dihedral annotation
# Add axis labels for dihedral annotation:
yy <- ymin + (ymax - ymin) * c(10.5/15, 11.5/15, 12.5/15)
axis(2, at=yy, labels=c("CCCC", "HCCC", "CCCH"), las=2, tck=.01, hadj=1.4)

dev.off()
