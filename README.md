
<!-- README.md is generated from README.Rmd. Please edit that file -->
PCMBase : Simulation and likelihood calculation of phylogenetic comparative methods
===================================================================================

Phylogenetic comparative models represent models of continuous trait data associated with the tips of a phylogenetic tree. Examples of such models are Gaussian continuous time branching stochastic processes such as Brownian motion (BM) and Ornstein-Uhlenbeck (OU), which regard the data at the tips of the tree as an observed (final) state of a Markov chain starting from an initial state at the root and evolving along the branches of the tree. The PCMBase R package provides a general framework for manipulating such models. This framework consists of an application programming interface for specifying data and model parameters, and efficient algorithms for simulating trait evolution under a model and calculating the likelihood of model parameters for an assumed model and trait data. The package implements a growing collection of models, which currently includes BM, OU, BM/OU with jumps, two-speed OU as well as mixed Gaussian models, in which different types of the above models can be associated with different branches of the tree. The PCMBase package is limited to trait-simulation and likelihood calculation of (mixed) Gaussian phylogenetic models. The PCMFit package provides functionality for ML and Bayesian fit of these models to tree and trait data.

Installation
------------

### CRAN

``` r
install.packages("PCMBase")
```

### Github

``` r
devtools::install_github("venelin/PCMBase")
```

Resources
---------

The user guides and technical reference for the library are available from the [PCMBase web-page](https://venelin.github.io/PCMBase/).

-   The [Getting started](https://venelin.github.io/PCMBase/articles/PCMBase.html)
-   The [The PCMBase parameterizations](https://venelin.github.io/PCMBase/articles/PCMParam.html) (in preparation)

The research article "Fast likelihood evaluation for multivariate phylogenetic comparative methods: the PCMBase R package" provides a general overview of PCMBase. The article is currently undergoing peer review for a publication and is available as a preprint from [arxiv](https://arxiv.org/abs/1809.09014).

The PCMBase source code is located in the [PCMBase github repository](https://github.com/venelin/PCMBase).

Feature requests, bugs, etc can be reported in the [PCMBase issues list](https://github.com/venelin/PCMBase/issues).
