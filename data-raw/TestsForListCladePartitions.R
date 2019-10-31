

  system.time({numAllPartsGt3 <- 0
  lstAllPartsGt3 <- list()
  K <- 0
  while(TRUE) {
    K <- K + 1L
    lstNewParts <- PCMTreeListCladePartitions(
      PCMBaseTestObjects$tree.ab, nNodes = K, minCladeSize = 3, countOnly = FALSE)
    lstAllPartsGt3 <- c(lstAllPartsGt3, lstNewParts)
    if(length(lstNewParts) == 0) {
      break
    }
  } })

  system.time(lstAllParts <- PCMTreeListAllPartitions(PCMBaseTestObjects$tree.ab, minCladeSize = 3))

  library(data.table)
  dtAllParts <- data.table(
    partition = c(list(PCMTreeNumTips(PCMBaseTestObjects$tree.ab) + 1L),
                  lstAllPartsGt3))
  dtAllParts[, K:=sapply(partition, length)]
  dtAllParts[, md5:=sapply(partition, function(p) {
    tree <- PCMBaseTestObjects$tree.ab
    PCMTreeSetPartition(tree, p)
    m <- PCMTreeMatrixNodesInSamePart(tree, nodes = seq_len(PCMTreeNumTips(tree)))
    digest::digest(m)
  })]

  dtAllPartsUnique <- dtAllParts[, list(K = K[1], partition = list(partition[[1]])), by = md5]


  numAllPartsGt3 <- 0
  lstAllPartsGt3 <- list()
  K <- 0
  while(TRUE) {
    K <- K + 1L
    numAllPartsGt3New <-
      PCMTreeListCladePartitions(PCMBaseTestObjects$tree_15_tips, nNodes = K, minCladeSize = 4, countOnly = TRUE)
    lstAllPartsGt3 <-
      c(lstAllPartsGt3,
        PCMTreeListCladePartitions(PCMBaseTestObjects$tree_15_tips, nNodes = K, minCladeSize = 4, countOnly = FALSE))
    if(numAllPartsGt3New == 0) {
      break
    } else {
      numAllPartsGt3 <- numAllPartsGt3 + numAllPartsGt3New
    }
  }

  lstAllParts <- PCMTreeListAllPartitions(PCMBaseTestObjects$tree_15_tips, minCladeSize = 4)




  compareLstParts <- function(tree, lstParts1, lstParts2) {
    for(i in seq_along(lstParts1)) {
      cat('.')
      PCMTreeSetPartition(tree, lstParts1[[i]])
      mat1 <- PCMTreeMatrixNodesInSamePart(tree)
      matchFound <- FALSE
      for(j in seq_along(lstParts2)) {
        PCMTreeSetPartition(tree, lstParts2[[j]])
        mat2 <- PCMTreeMatrixNodesInSamePart(tree)

        if(isTRUE(all.equal(mat1, mat2))) {
          matchFound <- TRUE
        }
      }
      if(!matchFound) {
        cat("\nNo match in lstParts2 for ", i)
      }
    }
  }

  compareLstParts(tree_15_tips2, lstAllPartsGt3, lstAllParts)
  compareLstParts(tree_15_tips2, lstAllParts, lstAllPartsGt3)

  PCMTreeSetPartition(tree_15_tips2, lstAllParts[[13]])
  PCMTreePlot(tree_15_tips2) + geom_tiplab() + geom_nodelab()

  browser()

  PCMTreeListCladePartitions(PCMBaseTestObjects$tree_15_tips, nNodes = 2, minCladeSize = 4, countOnly = FALSE)
