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
```R
library(devtools)
install_github("Melkiades/CampaRi")
library(CampaRi)
```

------------------------------------------------------------------------

### Manual ###

For more details, please refer to the main documentation of the original [CAMPARI documentation](http://campari.sourceforge.net/documentation.html).


### Usage - tests ###

All the tests can be found in the 'unusual' folder tests. The whole system must be still updated with the __testthat__ support and possibly a __vignette__ integration.


------------------------------------------------------------------------

## Citation

Please cite the SAPPHIRE algorithm paper: 

* Blöchliger, N., Vitalis, A. & Caflisch, A. A scalable algorithm to order and annotate continuous observations reveals the metastable states visited by dynamical systems. Comput. Phys. Commun. 184, 2446–2453 (2013).
<!---
[![doi](https://img.shields.io/badge/doi-10.1093/nar/gkv1507-green.svg?style=flat)](http://dx.doi.org/10.1093/nar/gkv1507) [![citation](https://img.shields.io/badge/cited%20by-18-green.svg?style=flat)](https://scholar.google.com.hk/scholar?oi=bibs&hl=en&cites=6029790855238928406) [![Altmetric](https://img.shields.io/badge/Altmetric-27-green.svg?style=flat)](https://www.altmetric.com/details/4919535)
-->

