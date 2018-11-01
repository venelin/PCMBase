.RunPCMBaseTests <- Sys.getenv("RunPCMBaseTests") == "yes"

if(.RunPCMBaseTests) {

  library(ape)
library(testthat)
library(PCMBase)

if(FALSE) {
  # MANUAL/VISUAL TEST ONLY
  set.seed(1)

  tree <- rtree(15)
  PCMTreeSetDefaultRegime(tree, 1)

  test_that("PCMTreeNodeTimes", identical(PCMTreeNodeTimes(tree), POUMM::nodeTimes(tree)))

  PCMTreeSetLabels(tree)
  par(mfrow=c(2, 3))
  plot(tree)
  nodelabels(text=tree$node.label)

  tree1 <- PCMTreeEvalNestedED(E(tree, 20))
  plot(tree1)
  nodelabels(text=tree1$node.label)

  tree2 <- PCMTreeEvalNestedED(E(E(tree, 20),28))
  plot(tree2)
  nodelabels(text=tree2$node.label)

  tree3 <- PCMTreeEvalNestedED(E(D(tree,27),29))
  plot(tree3)
  nodelabels(text=tree3$node.label)

  tree4 <- PCMTreeEvalNestedED(D(D(D(E(tree,22),17),26),25))
  plot(tree4)
  nodelabels(text = tree4$node.label)

  tree5 <- PCMTreeEvalNestedED(D(D(D(E(tree,20),22),10),27))
  plot(tree5)
  nodelabels(text = tree5$node.label)

  tree6 <- PCMTreeEvalNestedED(E(tree,1))
  plot(tree6)
  nodelabels(text = tree6$node.label)
}

if(FALSE) {
  set.seed(1)

  tree <- rtree(200)
  PCMTreeSetDefaultRegime(tree, 1)
  PCMTreeSetLabels(tree)

  PCMTreePlot(tree)

  preorderTree <- PCMTreePreorder(tree)
  tableAncestors <- PCMTreeTableAncestors(tree, preorder = preorderTree)

  # 0. Create a list of all possible clade-partitions of tree into clades not smaller than
  # minCladeSize[currentRecDepth]
  #
  # start from a list containing the trivial partition into one clade equal to the whole tree
  listCladePartitions <- list(integer(0))
  numPartNodes <- 1
  while(TRUE) {
    listNew <- PCMTreeListCladePartitions(tree, numPartNodes, 20, tableAncestors)
    if(length(listNew) == 0) {
      break
    } else {
      listCladePartitions <- c(listCladePartitions, listNew)
      numPartNodes <- numPartNodes + 1
    }
  }

  # 1. (fitsToClades) Perform a fit of each model-type to each clade
  # -> fitsToClades

  # 2. (fitsToTree)

  # 3. Choose best clade-partition
  bestPartition <- listCladePartitions[[746]]

  PCMTreeSetRegimes(tree, bestPartition)

  # 4. For each part in the best clade-partition, find its best sub-partition
  # with respect to the whole tree and the cut points and mapped regimes in the remaining
  # part.
  treePartRoot <- PCMTreeEvalNestedED(D(D(D(D(D(tree, 203), 223), 255), 285), 307))

  listCladePartitionsPart <- list(integer(0))
  numPartNodes <- 1
  while(TRUE) {
    listNew <- PCMTreeListCladePartitions(treePartRoot, numPartNodes, 12, tableAncestors[PCMTreeGetLabels(treePartRoot), PCMTreeGetLabels(treePartRoot)])
    if(length(listNew) == 0) {
      break
    } else {
      listCladePartitionsPart <- c(listCladePartitionsPart, listNew)
      numPartNodes <- numPartNodes + 1
    }
  }

  PCMTreePlot(tree)

  PCMTreeGetStartingNodesRegimes(tree)
  PCMTreeGetStartingNodesRegimes(PCMTreeSetRegimes(tree, c(PCMTreeGetStartingNodesRegimes(tree), PCMTreeMatchLabels(tree, PCMTreeGetLabels(treePartRoot)[listCladePartitionsPart[[15]]])), inplace = FALSE))

  c(PCMTreeGetStartingNodesRegimes(tree), PCMTreeMatchLabels(tree, PCMTreeGetLabels(treePartRoot)[listCladePartitionsPart[[2]]]))

  PCMTreePlot(PCMTreeSetRegimes(tree, c(PCMTreeGetStartingNodesRegimes(tree), PCMTreeMatchLabels(tree, PCMTreeGetLabels(treePartRoot)[listCladePartitionsPart[[15]]])), inplace = FALSE))

}
}
