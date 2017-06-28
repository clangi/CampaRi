<!---
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/BioinformaticsFMRP/TCGAbiolinks?branch=master&svg=true)](https://ci.appveyor.com/project/BioinformaticsFMRP/TCGAbiolinks)
[![codecov.io](https://codecov.io/github/BioinformaticsFMRP/TCGAbiolinks/coverage.svg?branch=master)](https://codecov.io/github/BioinformaticsFMRP/TCGAbiolinks?branch=master)
[![bioc](http://www.bioconductor.org/shields/downloads/TCGAbiolinks.svg)](http://bioconductor.org/packages/stats/bioc/TCGAbiolinks.html)
[![bioc](http://www.bioconductor.org/shields/years-in-bioc/TCGAbiolinks.svg)](http://bioconductor.org/packages/TCGAbiolinks/)
[![bioc](http://bioconductor.org/shields/availability/devel/TCGAbiolinks.svg)](http://bioconductor.org/packages/TCGAbiolinks/)

-->
[![Travis-CI Build Status](https://travis-ci.org/Melkiades/CampaRi.svg?branch=master)](https://travis-ci.org/Melkiades/CampaRi)

------------------------------------------------------------------------
# CampaRi - An R wrapper for fast trajectory analysis

Analysis algorithms extracted from the original 'campari' software package.
They consists in a kinetic annotation of the trajectory based on the minimum spanning tree constructed on the distances between snapshots. The fast algorithm is implemented on the basis of a modified version of the birch algorithm, while the slow one is based on a simple leader clustering. For more information please visit the [original website](http://campari.sourceforge.net/index.html)

### Installation ###

This code is available on this GitHub repository and therefore it can be downloaded directly and locally installed or, using the handy "devtools", it can be installed from R in the following simple way:

```R
library(devtools) # external libssl-dev is needed
install_github("Melkiades/CampaRi")

# Please note that you must install some libraries from Bioconductor (these are not automatically installed from dependencies)
# To do so: 
source("http://bioconductor.org/biocLite.R") 
biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")) 
install.packages("WGCNA") # this is the dependence that needs bioconductor

library(CampaRi)
```
Please note that to unlock some fundamental analysis function you need to install netcdf 4.1 or later library for backend fast data-handling.
This library must be installed --with-fortran bindings (it is possible to check this using the nc-config or nf-config commands (with --all option)).
Netcdf4 can be installed as following (this is the simple way):

```sh
# mac osx
brew update
brew install netcdf --with-fortran

# linux
sudo apt update -g
sudo apt install libnetcdf-dev --with-fortran
```

------------------------------------------------------------------------
### Usage - tests ###

We are currently writing how-to-use guides with examples and datasets. These vignettes are available internally in the package and can be browsed as following (please note that some external library is needed - e.g. pandoc).

```R
library(knitr)
library(rmarkdown)
install_github("Melkiades/CampaRi", build_vignettes= TRUE)

# TO BROWSE THEM after the installation
browseVignettes("CampaRi") 
```


------------------------------------------------------------------------
### Extensive Manual ###

For more details, please refer to the main documentation of the original [CAMPARI documentation](http://campari.sourceforge.net/documentation.html).






------------------------------------------------------------------------

## Citation

Please cite the SAPPHIRE algorithm paper: 

* Blöchliger, N., Vitalis, A. & Caflisch, A. A scalable algorithm to order and annotate continuous observations reveals the metastable states visited by dynamical systems. Comput. Phys. Commun. 184, 2446–2453 (2013).
<!---
[![doi](https://img.shields.io/badge/doi-10.1093/nar/gkv1507-green.svg?style=flat)](http://dx.doi.org/10.1093/nar/gkv1507) [![citation](https://img.shields.io/badge/cited%20by-18-green.svg?style=flat)](https://scholar.google.com.hk/scholar?oi=bibs&hl=en&cites=6029790855238928406) [![Altmetric](https://img.shields.io/badge/Altmetric-27-green.svg?style=flat)](https://www.altmetric.com/details/4919535)
-->

