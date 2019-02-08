---
title: "NEWS about the PCMBase R-package"
author: "Venelin Mitov"
date: "29 November, 2018"
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

sed -i "" 's/PCMTreeSetDefaultRegime/PCMTreeSetDefaultPartition/g' *
sed -i "" 's/PCMTreeGetStartingNodesRegimes/PCMTreeGetPartition/g' *
sed -i "" 's/PCMTreeSetRegimes/PCMTreeSetPartition/g' *
sed -i "" 's/PCMTreeUniqueRegimes/PCMTreeGetPartNames/g' *
sed -i "" 's/PCMTreeNumUniqueRegimes/PCMTreeNumParts/g' *
sed -i "" 's/PCMTreeGetRegimesForNodes/PCMTreeGetPartsForNodes/g' *
sed -i "" 's/PCMTreeGetTipsInRegime/PCMTreeGetTipsInPart/g' *
sed -i "" 's/PCMTreeMatchRegimesWithModel/PCMTreeGetRegimesForEdges/g' *
sed -i "" 's/PCMTreeDtNodeRegimes/PCMTreeDtNodeParts/g' *
sed -i "" 's/PCMTreeExtractBackboneRegimes/PCMTreeBackbonePartition/g' *
sed -i "" 's/PCMTreeMatrixNodesInSameRegime/PCMTreeMatrixNodesInSamePart/g' *

* 125 model types are made available within the package namespace. Hence, there
is no need to call PCMGenerateModelTypes() after loading the package. Call PCMModels()
after loading the package to see which are these model types. 

# PCMBase 1.2.8
* Added handling of standard error specified for each trait measurement. This is 
implemented as an additional matrix parameter called SE to the functions PCMLik, 
PCMSim, PCMVar, and PCMInfo.


# PCMBase 1.2.7
* First version of the package released on CRAN.

