library(ape)
library(testthat)
library(PCMBase)

if(FALSE) {
  # MANUAL/VISUAL TEST ONLY
  set.seed(1)

  tree <- rtree(15)
  PCMTreeSetDefaultRegime(tree, 1)

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
