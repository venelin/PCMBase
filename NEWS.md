---
title: "NEWS about the PCMBase R-package"
author: "Venelin Mitov"
date: "15 March, 2018"
output: html_document
---

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

