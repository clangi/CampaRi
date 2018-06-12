<!---
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/BioinformaticsFMRP/TCGAbiolinks?branch=master&svg=true)](https://ci.appveyor.com/project/BioinformaticsFMRP/TCGAbiolinks)
[![codecov.io](https://codecov.io/github/BioinformaticsFMRP/TCGAbiolinks/coverage.svg?branch=master)](https://codecov.io/github/BioinformaticsFMRP/TCGAbiolinks?branch=master)
[![bioc](http://www.bioconductor.org/shields/downloads/TCGAbiolinks.svg)](http://bioconductor.org/packages/stats/bioc/TCGAbiolinks.html)
[![bioc](http://www.bioconductor.org/shields/years-in-bioc/TCGAbiolinks.svg)](http://bioconductor.org/packages/TCGAbiolinks/)
[![bioc](http://bioconductor.org/shields/availability/devel/TCGAbiolinks.svg)](http://bioconductor.org/packages/TCGAbiolinks/)

[![codecov](https://codecov.io/gh/Melkiades/CampaRi/branch/master/graph/badge.svg)](https://codecov.io/gh/Melkiades/CampaRi)
-->
[![pipeline status](https://gitlab.com/CaflischLab/CampaRi/badges/master/pipeline.svg)](https://gitlab.com/CaflischLab/CampaRi/commits/master)
[![coverage report](https://gitlab.com/CaflischLab/CampaRi/badges/master/coverage.svg)](https://gitlab.com/CaflischLab/CampaRi/commits/master)
[![DOI](https://zenodo.org/badge/68593949.svg)](https://zenodo.org/badge/latestdoi/68593949)


------------------------------------------------------------------------
# CampaRi - An R wrapper for fast trajectory analysis

For any help with the package usage and installation please [chat with us](https://gitter.im/CampaRi_chatroom/).

Analysis algorithms extracted from the original 'campari' software package.
They consists in a kinetic annotation of the trajectory based on the minimum spanning tree constructed on the distances between snapshots. The fast algorithm is implemented on the basis of a modified version of the birch algorithm, while the slow one is based on a simple leader clustering. For more information please visit the [original website](http://campari.sourceforge.net/index.html)

### Installation ###

This code is available on this GitHub repository and therefore it can be downloaded directly and locally installed or, using the handy "devtools".
As we recently tested the folowing code on a fresh Ubuntu 16.04 we advise the following installations. 

Firstly from terminal we install the most modern versions of the following packages using apt:
```sh
# ubuntu 16.04 (remember to have R installed -> sudo apt install r-base*)
sudo apt update 
sudo apt install libnetcdff-dev # xenial: needed for light data handling
sudo apt install libarpack2-dev # xenial: needed for alternative to lapack spectral decomposition of matrices
sudo apt install libssl-dev # for devtools
sudo apt install libcurl4-gnutls-dev # for devtools
sudo apt install libtiff5-dev # for ClusterR
sudo apt install libgmp-dev # for ClusterR
# in the case of installation of original campari software: install_campari()
sudo apt install libfftw3-dev # for install_campari()

# some peculiarity for trusty as it had a different netcdf version and installation procedure
sudo apt-get install libnetcdff5 # trusty
sudo apt-get install libnetcdf-dev # trusty
sudo apt-get install libarpack2-dev # trusty

# osx # to-update 
brew update
brew install netcdf --with-fortran
brew install arpack
```

Secondly install extra-cran-packages using biocLite():
```R
source("http://bioconductor.org/biocLite.R") 
biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore", "PairViz")) 
```
Note: the installation of biocgenerics is ultimately failing on ubuntu 17.10 and only R > 3.4.3 can patch it up. We advise to stick with the LTS version 16.04 of ubuntu.

Finally install CampaRi package from github using devtools (note ggfortify just moved away from cran and must be installed from github). It is not needed anymore.

```R
library(devtools)
# install_github("sinhrks/ggfortify")  # no more needed

# Final step:
install_github("Melkiades/CampaRi")
```

Please note that to unlock some fundamental analysis functions you need to install netcdf 4.1 or later library for backend fast data-handling.
This library must be installed --with-fortran bindings (it is possible to check this using the nc-config or nf-config commands (with --all option)).

The arpack library is needed at this moment. It will be made optional in the future.

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

* Bloechliger, N., Vitalis, A. & Caflisch, A. A scalable algorithm to order and annotate continuous observations reveals the metastable states visited by dynamical systems. Comput. Phys. Commun. 184, 2446–2453 (2013).
<!---
[![doi](https://img.shields.io/badge/doi-10.1093/nar/gkv1507-green.svg?style=flat)](http://dx.doi.org/10.1093/nar/gkv1507) [![citation](https://img.shields.io/badge/cited%20by-18-green.svg?style=flat)](https://scholar.google.com.hk/scholar?oi=bibs&hl=en&cites=6029790855238928406) [![Altmetric](https://img.shields.io/badge/Altmetric-27-green.svg?style=flat)](https://www.altmetric.com/details/4919535)
-->

