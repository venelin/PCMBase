
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Travis build status](https://travis-ci.org/venelin/PCMBase.svg?branch=master)](https://travis-ci.org/venelin/PCMBase) [![Coverage status](https://codecov.io/gh/venelin/PCMBase/branch/master/graph/badge.svg)](https://codecov.io/github/venelin/PCMBase?branch=master) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/PCMBase?color=blue)](https://cran.r-project.org/package=PCMBase) [![Downloads](http://cranlogs.r-pkg.org/badges/PCMBase?color=blue)](https://cran.r-project.org/package=PCMBase) [![Latest release DOI](https://zenodo.org/badge/144310133.svg)](https://zenodo.org/badge/latestdoi/144310133)

PCMBase : Simulation and likelihood calculation of phylogenetic comparative methods
===================================================================================

Phylogenetic comparative methods represent models of continuous trait data associated with the tips of a phylogenetic tree. Examples of such models are Gaussian continuous time branching stochastic processes such as Brownian motion (BM) and Ornstein-Uhlenbeck (OU) processes, which regard the data at the tips of the tree as an observed (final) state of a Markov process starting from an initial state at the root and evolving along the branches of the tree. The PCMBase R package provides a general framework for manipulating such models. This framework consists of an application programming interface for specifying data and model parameters, and efficient algorithms for simulating trait evolution under a model and calculating the likelihood of model parameters for an assumed model and trait data. The package implements a growing collection of models, which currently includes BM, OU, BM/OU with jumps, two-speed OU as well as mixed Gaussian models, in which different types of the above models can be associated with different branches of the tree. The PCMBase package is limited to trait-simulation and likelihood calculation of (mixed) Gaussian phylogenetic models. The PCMFit package provides functionality for ML and Bayesian fit of these models to tree and trait data.

Installation
------------

### An optional but highly recommended dependency

The function `PCMTreePlot` in the package is implemented based on the R-package ggtree, which is not on CRAN. It is highly recommended to install this package in order to be able to visualize trees with colored parts corresponding to defferent evolutionary regimes. If ggtree is not installed, the package will fail to run examples and generate the vignettes. At the time of writing this documentation, ggtree can be installed from bioconductor through the following code (if that does not work, check the [ggtree home page](https://guangchuangyu.github.io/software/ggtree/)):

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ggtree", version = "3.8")
```

### Installing PCMBase from CRAN

A stable but possibly old version of PCMBase is available on CRAN and can be installed with this command:

``` r
install.packages("PCMBase")
```

### Github

The newest but possibly less stable and tested version of the package can be installed using:

``` r
devtools::install_github("venelin/PCMBase")
```

Resources
---------

The user guides and technical reference for the library are available from the [PCMBase web-page](https://venelin.github.io/PCMBase/).

-   The [Getting started](https://venelin.github.io/PCMBase/articles/PCMBase.html)
-   The [The PCMBase parameterizations](https://venelin.github.io/PCMBase/articles/PCMParam.html) (in preparation)
-   The [Creating a Custom Model in the PCMBase Framework](https://venelin.github.io/PCMBase/articles/PCMCreateModel.html)

The research article "Fast likelihood evaluation for multivariate phylogenetic comparative methods: the PCMBase R package" provides a general overview of PCMBase. The article is currently undergoing peer review for a publication and is available as a preprint from [arxiv](https://arxiv.org/abs/1809.09014).

The PCMBase source code is located in the [PCMBase github repository](https://github.com/venelin/PCMBase).

Feature requests, bugs, etc can be reported in the [PCMBase issues list](https://github.com/venelin/PCMBase/issues).

Citing PCMBase
==============

To give credit to the PCMBase package in a publication, please cite the following article:

Mitov, V., Bartoszek, K., Asimomitis, G., & Stadler, T. (2018, September 24). Fast likelihood evaluation for multivariate phylogenetic comparative methods: the PCMBase R package. arXiv.org. <https://arxiv.org/abs/1809.09014>.

Used R-packages
===============

The PCMBase R-package uses the following 3rd party R-packages:

-   For tree processing in R: ape v5.3 (Paradis et al. 2018), data.table v1.12.0 (Dowle and Srinivasan 2019);
-   For algebraic manipulation: expm v0.999.4 (Goulet et al. 2018), mvtnorm v1.0.10 (Genz et al. 2018);
-   For plotting: ggtree v1.14.6 (Yu and Lam 2019), ggplot2 v3.1.0 (Wickham et al. 2018);
-   For unit-testing: testthat v2.0.1 (Wickham 2018), covr v3.2.1 (Hester 2018);
-   For documentation and web-site generation: roxygen2 v6.1.1 (Wickham, Danenberg, and Eugster 2018), pkgdown v1.3.0 (Wickham and Hesselberth 2018);

Licence and copyright
=====================

Copyright 2016-2019 Venelin Mitov

Source code to PCMBase is made available under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. PCMBase is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

References
==========

Dowle, Matt, and Arun Srinivasan. 2019. *Data.table: Extension of ‘Data.frame‘*. <https://CRAN.R-project.org/package=data.table>.

Genz, Alan, Frank Bretz, Tetsuhisa Miwa, Xuefei Mi, and Torsten Hothorn. 2018. *Mvtnorm: Multivariate Normal and T Distributions*. <https://CRAN.R-project.org/package=mvtnorm>.

Goulet, Vincent, Christophe Dutang, Martin Maechler, David Firth, Marina Shapira, and Michael Stadelmann. 2018. *Expm: Matrix Exponential, Log, ’Etc’*. <https://CRAN.R-project.org/package=expm>.

Hester, Jim. 2018. *Covr: Test Coverage for Packages*. <https://CRAN.R-project.org/package=covr>.

Paradis, Emmanuel, Simon Blomberg, Ben Bolker, Joseph Brown, Julien Claude, Hoa Sien Cuong, Richard Desper, et al. 2018. *Ape: Analyses of Phylogenetics and Evolution*. <https://CRAN.R-project.org/package=ape>.

Wickham, Hadley. 2018. *Testthat: Unit Testing for R*. <https://CRAN.R-project.org/package=testthat>.

Wickham, Hadley, and Jay Hesselberth. 2018. *Pkgdown: Make Static Html Documentation for a Package*. <https://CRAN.R-project.org/package=pkgdown>.

Wickham, Hadley, Winston Chang, Lionel Henry, Thomas Lin Pedersen, Kohske Takahashi, Claus Wilke, and Kara Woo. 2018. *Ggplot2: Create Elegant Data Visualisations Using the Grammar of Graphics*. <https://CRAN.R-project.org/package=ggplot2>.

Wickham, Hadley, Peter Danenberg, and Manuel Eugster. 2018. *Roxygen2: In-Line Documentation for R*. <https://CRAN.R-project.org/package=roxygen2>.

Yu, Guangchuang, and Tommy Tsan-Yuk Lam. 2019. *Ggtree: An R Package for Visualization and Annotation of Phylogenetic Trees with Their Covariates and Other Associated Data*. <https://guangchuangyu.github.io/software/ggtree>.
