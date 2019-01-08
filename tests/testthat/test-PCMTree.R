library(testthat)
context("PCMTree")

library(PCMBase)

if(PCMBaseIsADevRelease(numVersionComponents = 3)) {

  load("testobjects.RData")

  library(ape)


  set.seed(1)
  # number of tips
  N <- 40

  # tree with one regime
  tree.a <- rtree(N)
  PCMTreeSetDefaultRegime(tree.a, model.a.1)
  PCMTreeSetLabels(tree.a)
  #PCMTreePlot(tree.a)

  test_that(
    "Tree with 40 tips and one regime", {
      expect_identical(PCMTreeGetParent(tree.a, 71L), 52L)
      expect_identical(PCMTreeGetLabels(tree.a)[52L], "52")
      expect_identical(PCMTreeGetRegimesForNodes(tree.a, "71"), "a")
      expect_identical(PCMTreeNumUniqueRegimes(tree.a), 1L)
      expect_identical(PCMTreeGetStartingNodesRegimes(tree.a), c(a = 41L))
    })

  tree.ab <- tree.a
  PCMTreeSetRegimes(tree.ab, nodes = N + 31, regimes = c("a", "b"))
  #PCMTreePlot(tree.ab)

  test_that(
    "Tree with 40 tips and two regimes", {
      expect_identical(PCMTreeGetParent(tree.ab, 71L), 52L)
      expect_identical(PCMTreeGetLabels(tree.ab)[52L], "52")
      expect_identical(PCMTreeGetRegimesForNodes(tree.ab, "71"), "b")
      expect_identical(PCMTreeGetStartingNodesRegimes(tree.ab), c(a = 41L, b = 71L))
    })


  set.seed(1)

  tree <- rtree(15)
  PCMTreeSetDefaultRegime(tree, 1)

  test_that("PCMTreeNodeTimes", {
    expect_equivalent(
      PCMTreeNodeTimes(tree),
      c(0.666725094430149, 1.03572271834128, 1.47726051113568, 1.61899162991904, 2.10669805249199, 3.57368590473197, 3.26124938833527, 3.56297762691975, 2.59485166589729, 1.87603989476338, 2.49500915128738, 2.94890801399015, 2.41463545593433, 2.32295331568457, 2.70654317038134, 0, 0.267220668727532, 0.653334761271253, 1.136911514448, 0.599565825425088, 1.09310713247396, 1.27932473388501, 1.94779147207737, 2.74203133280389, 2.84997495869175, 1.15260213706642, 1.68232171726413, 2.47167794895358, 1.84533369354904))
  })

  PCMTreeSetLabels(tree)

  # par(mfrow=c(2, 3))
  # plot(tree); nodelabels(text=tree$node.label)

  tree1 <- PCMTreeEvalNestedEDxOnTree("E(tree, 20)", tree)
  # plot(tree1); nodelabels(text=tree1$node.label)
  test_that("Test PCMTreeEvalNestedEDxOnTree 1: E(tree, 20)", {
    expect_identical(PCMTreeGetLabels(tree1)[PCMTreeNumTips(tree1) + 1], "20")
  })

  tree2 <- PCMTreeEvalNestedEDxOnTree("E(E(tree, 20),28)", tree)
  # plot(tree2); nodelabels(text=tree2$node.label)
  test_that("Test PCMTreeEvalNestedEDxOnTree 2: E(E(tree, 20),28)", {
    expect_identical(PCMTreeGetLabels(tree2)[PCMTreeNumTips(tree2) + 1], "28")
  })

  tree3 <- PCMTreeEvalNestedEDxOnTree("E(D(tree,27),29)", tree)
  #plot(tree3); nodelabels(text=tree3$node.label)
  test_that("Test PCMTreeEvalNestedEDxOnTree 3: E(D(tree,27),29)", {
    expect_identical(PCMTreeGetLabels(tree3)[PCMTreeNumTips(tree3) + 1], "29")
  })

  tree4 <- PCMTreeEvalNestedEDxOnTree("D(D(D(E(tree,22),17),26),25)", tree)
  # plot(tree4); nodelabels(text = tree4$node.label)
  test_that("Test PCMTreeEvalNestedEDxOnTree 4: D(D(D(E(tree,22),17),26),25)", {
    expect_identical(PCMTreeGetLabels(tree4)[PCMTreeNumTips(tree4) + 1], "22")
    expect_error(PCMTreeGetParent(tree4, "8"))
    expect_identical(PCMTreeGetParent(tree4, PCMTreeMatchLabels(tree4, "8")), PCMTreeMatchLabels(tree4, "24"))
    expect_identical(PCMTreeGetParent(tree4, PCMTreeMatchLabels(tree4, "24")), PCMTreeMatchLabels(tree4, "23"))
  })

  #tree5 <- PCMTreeEvalNestedEDxOnTree("D(D(D(E(tree,20),22),10),27)", tree)
  #plot(tree5); nodelabels(text = tree5$node.label)

  #tree6 <- PCMTreeEvalNestedEDxOnTree("D(tree,1)", tree)
  #plot(tree6); nodelabels(text = tree6$node.label)


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
