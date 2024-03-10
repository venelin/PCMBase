
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Coverage
status](https://codecov.io/gh/venelin/PCMBase/branch/master/graph/badge.svg)](https://app.codecov.io/github/venelin/PCMBase?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/PCMBase?color=blue)](https://cran.r-project.org/package=PCMBase)
[![Downloads](http://cranlogs.r-pkg.org/badges/PCMBase?color=blue)](https://cran.r-project.org/package=PCMBase)

# PCMBase : Simulation and likelihood calculation of phylogenetic comparative methods

Phylogenetic comparative methods represent models of continuous trait
data associated with the tips of a phylogenetic tree. Examples of such
models are Gaussian continuous time branching stochastic processes such
as Brownian motion (BM) and Ornstein-Uhlenbeck (OU) processes, which
regard the data at the tips of the tree as an observed (final) state of
a Markov process starting from an initial state at the root and evolving
along the branches of the tree. The PCMBase R package provides a general
framework for manipulating such models. This framework consists of an
application programming interface for specifying data and model
parameters, and efficient algorithms for simulating trait evolution
under a model and calculating the likelihood of model parameters for an
assumed model and trait data. The package implements a growing
collection of models, which currently includes BM, OU, BM/OU with jumps,
two-speed OU as well as mixed Gaussian models, in which different types
of the above models can be associated with different branches of the
tree. Note that the PCMBase package does not implement model inference.
Due to the enormous variety of models and possible model inference
methods, this functionality is delegated to other packages that, taking
advantage of PCMBase’s fast likelihood calculation, can implement
maximum likelihood (ML) or Bayesian inference methods. For example, the
[PCMFit package](https://venelin.github.io/PCMFit/) provides
heuristic-based ML fit of (mixed) Gaussian phylogenetic models (Mitov,
Bartoszek, and Stadler 2019).

## Installation

### An optional but highly recommended dependency

The function `PCMTreePlot` in the package is implemented based on the
R-package ggtree, which is not on CRAN. It is highly recommended to
install this package in order to visualize trees with colored parts
corresponding to different evolutionary regimes. If ggtree is not
installed, the package will fail to run examples and generate the
vignettes. At the time of writing this documentation, ggtree can be
installed from bioconductor through the following code (if that does not
work, check the [ggtree home
page](https://guangchuangyu.github.io/software/ggtree/)):

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ggtree", version = "3.8")
```

### Installing PCMBase from CRAN

A stable but possibly old version of PCMBase is available on CRAN and
can be installed with this command:

``` r
install.packages("PCMBase")
```

### Github

The newest but possibly less stable and tested version of the package
can be installed using:

``` r
devtools::install_github("venelin/PCMBase")
```

## Resources

The user guides and technical reference for the library are available
from the [PCMBase web-page](https://venelin.github.io/PCMBase/).

-   The [Getting
    started](https://venelin.github.io/PCMBase/articles/PCMBase.html)
    guide
-   The [PCMBase
    parameterizations](https://venelin.github.io/PCMBase/articles/PCMParam.html)
    guide
-   The [Tracing the likelihood calculation of a Gaussian
    mode](https://venelin.github.io/PCMBase/articles/PCMTracePruning.html)
    guide
-   The [Creating a Custom Model in the PCMBase
    Framework](https://venelin.github.io/PCMBase/articles/PCMCreateModel.html)
    guide

The research article “Fast likelihood calculation for multivariate
Gaussian phylogenetic models with shifts”, published in *Theoretical
Population Biology* provides a thorough description of the likelihood
calculation algorithm currently implemented in PCMBase. Appendix A of
this article gives an overview of the modular structure and the features
of the package.

The PCMBase source code is located in the [PCMBase github
repository](https://github.com/venelin/PCMBase).

Feature requests, bugs, etc can be reported in the [PCMBase issues
list](https://github.com/venelin/PCMBase/issues).

# Citing PCMBase

To give credit to the PCMBase package in a publication, please cite the
following article:

Mitov, V., Bartoszek, K., Asimomitis, G., & Stadler, T. (2019). Fast
likelihood calculation for multivariate Gaussian phylogenetic models
with shifts. Theor. Popul. Biol.
<https://doi.org/10.1016/j.tpb.2019.11.005>

# Used R-packages

The PCMBase R-package uses the following 3rd party R-packages:

-   For tree processing in R: ape v5.5 (Paradis et al. 2018), data.table
    v1.14.2 (Dowle and Srinivasan 2019);
-   For algebraic manipulation: expm v0.999.6 (Goulet et al. 2018),
    mvtnorm v1.1.3 (Genz et al. 2018);
-   For plotting: ggtree v3.2.1 (Yu and Lam 2019), ggplot2 v3.3.5
    (Wickham et al. 2018);
-   For unit-testing: testthat v3.1.1 (Wickham 2018), covr v3.5.1
    (Hester 2018);
-   For documentation and web-site generation: roxygen2 v7.1.2 (Wickham,
    Danenberg, and Eugster 2018), pkgdown v2.0.1 (Wickham and
    Hesselberth 2018);

# Licence and copyright

Copyright 2016-2021 Venelin Mitov

Source code to PCMBase is made available under the terms of the GNU
General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version.
PCMBase is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
more details.

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-R-data.table" class="csl-entry">

Dowle, Matt, and Arun Srinivasan. 2019. *Data.table: Extension of
‘Data.frame‘*. <https://CRAN.R-project.org/package=data.table>.

</div>

<div id="ref-R-mvtnorm" class="csl-entry">

Genz, Alan, Frank Bretz, Tetsuhisa Miwa, Xuefei Mi, and Torsten Hothorn.
2018. *Mvtnorm: Multivariate Normal and t Distributions*.
<https://CRAN.R-project.org/package=mvtnorm>.

</div>

<div id="ref-R-expm" class="csl-entry">

Goulet, Vincent, Christophe Dutang, Martin Maechler, David Firth, Marina
Shapira, and Michael Stadelmann. 2018. *Expm: Matrix Exponential, Log,
’Etc’*. <https://CRAN.R-project.org/package=expm>.

</div>

<div id="ref-R-covr" class="csl-entry">

Hester, Jim. 2018. *Covr: Test Coverage for Packages*.
<https://CRAN.R-project.org/package=covr>.

</div>

<div id="ref-Mitov:2019ci" class="csl-entry">

Mitov, Venelin, Krzysztof Bartoszek, and Tanja Stadler. 2019. “<span
class="nocase">Automatic generation of evolutionary hypotheses using
mixed Gaussian phylogenetic models</span>.” *Proceedings of the National
Academy of Sciences of the United States of America* 35 (August):
201813823.

</div>

<div id="ref-R-ape" class="csl-entry">

Paradis, Emmanuel, Simon Blomberg, Ben Bolker, Joseph Brown, Julien
Claude, Hoa Sien Cuong, Richard Desper, et al. 2018. *Ape: Analyses of
Phylogenetics and Evolution*. <https://CRAN.R-project.org/package=ape>.

</div>

<div id="ref-R-testthat" class="csl-entry">

Wickham, Hadley. 2018. *Testthat: Unit Testing for r*.
<https://CRAN.R-project.org/package=testthat>.

</div>

<div id="ref-R-ggplot2" class="csl-entry">

Wickham, Hadley, Winston Chang, Lionel Henry, Thomas Lin Pedersen,
Kohske Takahashi, Claus Wilke, and Kara Woo. 2018. *Ggplot2: Create
Elegant Data Visualisations Using the Grammar of Graphics*.
<https://CRAN.R-project.org/package=ggplot2>.

</div>

<div id="ref-R-roxygen2" class="csl-entry">

Wickham, Hadley, Peter Danenberg, and Manuel Eugster. 2018. *Roxygen2:
In-Line Documentation for r*.
<https://CRAN.R-project.org/package=roxygen2>.

</div>

<div id="ref-R-pkgdown" class="csl-entry">

Wickham, Hadley, and Jay Hesselberth. 2018. *Pkgdown: Make Static HTML
Documentation for a Package*.
<https://CRAN.R-project.org/package=pkgdown>.

</div>

<div id="ref-R-ggtree" class="csl-entry">

Yu, Guangchuang, and Tommy Tsan-Yuk Lam. 2019. *Ggtree: An r Package for
Visualization and Annotation of Phylogenetic Trees with Their Covariates
and Other Associated Data*.
<https://guangchuangyu.github.io/software/ggtree/>.

</div>

</div>
