sed -i "" 's/PCMTreeSetDefaultRegime/PCMTreeSetDefaultPartition/g' *
sed -i "" 's/PCMTreeGetStartingNodesRegimes/PCMTreeGetPartition/g' *
# sed -i ""'s/StartingNodesRegimes/Partition/g' *
sed -i "" 's/PCMTreeSetRegimes/PCMTreeSetPartition/g' *
sed -i "" 's/PCMTreeUniqueRegimes/PCMTreeGetPartNames/g' *
sed -i "" 's/PCMTreeNumUniqueRegimes/PCMTreeNumParts/g' *
sed -i "" 's/PCMTreeGetRegimesForNodes/PCMTreeGetPartsForNodes/g' *
sed -i "" 's/PCMTreeGetTipsInRegime/PCMTreeGetTipsInPart/g' *
sed -i "" 's/PCMTreeMatchRegimesWithModel/PCMTreeGetRegimesForEdges/g' *
sed -i "" 's/PCMTreeDtNodeRegimes/PCMTreeDtNodeParts/g' *
sed -i "" 's/PCMTreeExtractBackboneRegimes/PCMTreeBackbonePartition/g' *
sed -i "" 's/PCMTreeMatrixNodesInSameRegime/PCMTreeMatrixNodesInSamePart/g' *
 
