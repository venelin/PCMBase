---
title: "NEWS about the PCMBase R-package"
author: "Venelin Mitov"
date: "15 November, 2019"
output: html_document
---

# PCMBase 1.2.11
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

