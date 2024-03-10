---
title: "NEWS about the PCMBase R-package"
author: "Venelin Mitov"
date: "17 Nov, 2021"
output: html_document
---

# PCMBase 1.2.14
- Fixed CRAN check issues in inst/CITATION (replacing calls to personList with c; replacing citEntry with bibentry).
- Fixed CRAN check issues regarding abind code in vignettes signaled by Prof. Brian Ripley.

# PCMBase 1.2.13
- Fixed CRAN check issues regarding ggtree code in vignettes signaled by Prof. Brian Ripley.

# PCMBase 1.2.12 
- Added function PCMTreeInsertTipsOrSingletons.
- Added possibility to specify metaI as a character string (function name) rather.
So far, it was only possible to specify this as a function object or as the list
resulting from calling the function on the tree, data and PCM model.
- Added an operator + for PCMTree objects. This is analogical to the operator +
for phylo objects but maintains the edge.part and part.regime members 
(experimental). 
- Removed argument parametrizations from PCMGenerateModelTypes. Now it's possible
to specify the type of parametrization (all/default) separately for each model 
type.
- Added automatic generation of White noise model types. 
- Improved man pages.

# PCMBase 1.2.11
- Updated package citation pointing to the recently published article: Mitov et al. 2019. Fast likelihood calculation for multivariate Gaussian phylogenetic models with shifts. Theoretical Population Biology. https://doi.org/10.1016/j.tpb.2019.11.005.
- Fixed an issue causing a CRAN-check failure on one of the flavors: https://www.r-project.org/nosvn/R.check/r-devel-linux-x86_64-debian-clang/PCMBase-00check.html
- Added functions PCMListMembers, MatchListMembers, PCMSetAttribute, PCMGetAttribute.
- Added examples and improved documentation of the functions PCMLik and UpperTriFactor. 
- Fixed a bug in the parsing of error messages (attribute error of the result
from calling PCMLik).
- Fixed a bug in PCMTreeListCladePartitions. This bug was manifesting only when nNodes > 1. Hence, previous usages in PCMFit R-package are not affected.
- Added argument countOnly to PCMTreeListCladePartitions.
- Added new functions PCMListMembers, MatchListMembers, PCMGetAttribute, 
PCMSetAttribute.
- Minor enhancement in the MixedGaussian constructor allowing to provide
a list of PCM objects as modelTypes argument.
- Added an argument diagOnly to PCMVar. 
- Added an argument enclos in PCMParamLocateInShortVector.
- Additions to the documentation.

# PCMBase 1.2.10
- Speedup of PCMApplyTransformation in the case of mixed Gaussian models. 
- Added functions PCMExtractDimensions and PCMExtractRegimes
- Added two new vignettes: PCMCreateModel.Rmd and PCMTracePruning.Rmd.
- Added a new model 'BMdrift'.
- Finished parametrization vignette. 
- Added functions for printing a PCM object as a table - see corresponding 
section in the Getting started guide.
- Added functions for tracing a likelihood calculation - see ?PCMLikTrace.
- Some new runtime options were added - see ?PCMOptions.
- Updated documentation pages.
- Several bugs were found and fixed.


# PCMBase 1.2.9
* Important changes in the PCMTree module: 
- Renamed member 'edge.regime' in the phylo object to 'edge.part'. The word 
'regime' is reserved uniquely to identify a regime of an evolutionary model. An
intuitive synonyme for 'regime' is 'color'. The word 'part' stays for a 
monophyletic or a paraphyletic group of nodes (and the branches ending at these 
nodes) in a tree. It is possible to have several parts that evolve under the 
same model regime and, terefore, should be shown with with the same color.

* new model types are made available within the package namespace. Hence, there
is no need to call PCMGenerateModelTypes() after loading the package. Call PCMModels()
after loading the package to see which are these model types. 

* Updated tests and examples due to a change in the default random generator type (see https://developer.r-project.org/blosxom.cgi/R-devel/NEWS/2019/02/28#n2019-02-28).

# PCMBase 1.2.8
* Added handling of standard error specified for each trait measurement. This is 
implemented as an additional matrix parameter called SE to the functions PCMLik, 
PCMSim, PCMVar, and PCMInfo.


# PCMBase 1.2.7
* First version of the package released on CRAN.

