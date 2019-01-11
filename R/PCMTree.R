# Copyright 2018 Venelin Mitov
#
# This file is part of PCMBase.
#
# PCMBase is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PCMBase is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PCMBase.  If not, see <http://www.gnu.org/licenses/>.

#' Wrapper for length(tree$tip.label)
#' @param tree a phylo object
#' @return the number of tips in tree
#' @export
PCMTreeNumTips <- function(tree) {
  length(tree$tip.label)
}

#' Number of all nodes in a tree
#'
#' @details Wrapper for nrow(tree$edge) + 1
#' @param tree a phylo object
#' @return the number of nodes in tree including root, internal and tips.
#' @export
PCMTreeNumNodes <- function(tree) {
  nrow(tree$edge) + 1L
}

#' Set a default edge.regime member ot the passed tree object
#' @param tree a phylo object
#' @param regime a character, an integer or PCM model object
#'
#' @description This function sets or overwrites the current member edge.regime
#' in tree with one of the following:
#' \itemize{
#' \item{rep(regime[1], length(tree$edge.length))}{if regime is a character or an integer}
#' \item{rep(PCMRegimes(regime)[1], length(tree$edge.length))}{if regime is a PCM model}
#' } Note that the function modifies the passed tree object inplace.
#' @return This function does not return a value but has a side effect on the passed
#' tree object.
#' @export
PCMTreeSetDefaultRegime <- function(tree, regime) {
  if(is.PCM(regime)) {
    regime <- PCMRegimes(regime)
  }
  eval(substitute(tree$edge.regime <- rep(regime[1], length(tree$edge.length))), parent.frame())
}

#' Assign regimes on a tree given a set of starting branches
#'
#' @details It is assumed that each regime "paints" a linked subset of branches
#' on a tree. Thus, each regime is fully described by its starting branch. The
#' descendant branches inherit this regime until reaching a tip or a node that is
#' present in the \code{nodes} parameter.
#'
#' @param tree a phylo object
#' @param nodes a character vector containing tip or node labels or an integer
#' vector denoting tip or internal nodes in tree - the regimes change at the
#' start of the branches leading to these nodes.
#' @param regimes an integer or character vector of length equal to
#' length(nodes) + 1 containing the regime-names to be assigned for each regime.
#' If NULL the regime names will be the integers 1:(length(nodes) + 1).
#' @param inplace a logical indicating if the change should be done to the tree
#' in the calling environment (TRUE) or a copy of the tree with set edge.regime
#' member should be returned (FALSE). Default is TRUE.
#' @return If inplace is TRUE nothing, otherwise a copy of the tree with set
#' edge.regime member.
#'
#' @seealso \code{\link{PCMTreeGetStartingNodesRegimes}}
#' @export
PCMTreeSetRegimes <- function(tree, nodes, regimes = as.integer(1:(length(nodes) + 1)), inplace=TRUE) {
  if(!inherits(tree, "phylo")) {
    stop("ERR:026d0:PCMBase:PCMTree.R:PCMTreeSetRegimes:: argument tree should be a phylo.")
  }
  if(!is.null(regimes)) {
    if(length(regimes) != length(nodes) + 1L ||
       length(unique(regimes)) < length(regimes)) {
      stop("ERR:026d1:PCMBase:PCMTree.R:PCMTreeSetRegimes:: regimes should be a character or integer vector of length equal to length(nodes) + 1 and not have duplicated elements.")
    }
  }
  if(is.character(nodes)) {
    nodes <- match(nodes, PCMTreeGetLabels(tree))
    if(any(is.na(nodes))) {
      stop("ERR:026d2:PCMBase:PCMTree.R:PCMTreeSetRegimes:: if nodes is a character vector it should be a subset of PCMTreeGetLabels(tree), but some of the elements in nodes could not be matched.")
    }
  }
  preorder <- PCMTreePreorder(tree)
  edge.regime <- rep(1, nrow(tree$edge))
  nextRegime <- 2
  N <- PCMTreeNumTips(tree)

  endingAt <- order(rbind(tree$edge, c(0, N+1))[, 2])


  for(ei in preorder) {
    i <- tree$edge[ei, 2]
    j <- tree$edge[ei, 1]
    ej <- endingAt[j]

    if(i %in% nodes) {
      edge.regime[ei] <- nextRegime
      nextRegime <- nextRegime + 1
    } else if(j != N+1) {
      edge.regime[ei] <- edge.regime[ej]
    }
  }

  if(!is.null(regimes)) {
    edge.regime <- regimes[edge.regime]
  }
  if(inplace) {
    eval(substitute(tree[["edge.regime"]] <- edge.regime), parent.frame())
  } else {
    tree$edge.regime <- edge.regime
    tree
  }
}

#' Set tip and internal node labels in a tree
#' @param tree a phylo object
#' @param labels a character vector in the order 1:PCMTreeNumNodes(tree) as denoted in the
#' tree$edge matrix.
#' @param inplace a logical indicating if the change should be done in place on the object
#' in the calling environment (in this case tree must not be a temporary object, e.g.
#' returned by another function call). Default is TRUE.
#' @return if inplace = FALSE, a copy of tree with set or modified tree$tip.label and tree$node.label, otherwise nothing.
#' @export
PCMTreeSetLabels <- function(tree, labels = as.character(1:PCMTreeNumNodes(tree)), inplace = TRUE) {
  N <- PCMTreeNumTips(tree)
  M <- PCMTreeNumNodes(tree)
  labels2 <- labels
  if(inplace) {
    eval(substitute({
      tree$tip.label <- labels2[1:N]
      tree$node.label <- labels2[(N + 1):M]
    }), parent.frame())
  } else {
    tree$tip.label <- labels2[1:N]
    tree$node.label <- labels2[(N + 1):M]
    tree
  }
}

#' Get a vector of the tip and node labels in a tree
#' @param tree a phylo object
#' @return a character vector
#' @export
PCMTreeGetLabels <- function(tree) {
  if(!inherits(tree, "phylo")) {
    stop("ERR:02650:PCMBase:PCMTree.R:PCMTreeGetLabels:: argument tree should be a phylo.")
  }
  if(is.null(tree$node.label)) {
    stop("ERR:02651:PCMBase:PCMTree.R:PCMTreeGetLabels:: the tree has no node.label member assigned.")
  }
  c(tree$tip.label, tree$node.label)
}

#' Get the node numbers associated with tip- or node-labels in a tree
#' @param tree a phylo object
#' @param labels a character vector with valid tip or node labels from tree
#' @return an integer vector giving the tip- or node- integer indices corresponding to labels.
#' @export
PCMTreeMatchLabels <- function(tree, labels) {
  allLabels <- PCMTreeGetLabels(tree)
  m <- match(labels, allLabels)
  if(any(is.na(m))) {
    stop("ERR:02660:PCMBase:PCMTree.R:PCMTreeMatchLabels:: some of the labels not found in PCMTreeGetLabels(tree).")
  }
  m
}

#' Get the starting branch' nodes for each regime on a tree
#'
#' @details We call a starting branch the first branch from the root to the tips
#' with a given regime. A starting node is the node at which a starting branch
#' ends.
#' @param tree a phylo object with an edge.regime member denoting regimes. The
#' function assumes that each regime covers a linked set of branches on the tree.
#' @param preorder an integer vector of row-indices in tree$edge as returned by
#' \code{\link{PCMTreePreorder}}. Defaults to \code{PCMTreePreorder(tree)}.
#' Specifying this argument may improve performance if PCMTreePreorder had to be
#' called earlier.
#' @return an integer with elements equal to the starting nodes for each regime in
#' \code{regimes}.
#' @seealso \code{\link{PCMTreeSetRegimes}}
#' @export
PCMTreeGetStartingNodesRegimes <- function(tree, preorder = PCMTreePreorder(tree)) {
  if(!inherits(tree, "phylo")) {
    stop("ERR:026g0:PCMBase:PCMTree.R:PCMTreeGetStartNodesRegimes:: argument tree should be a phylo.")
  }

  regimes <-  PCMTreeUniqueRegimes(tree)

  N <- PCMTreeNumTips(tree)
  nodes <- integer(length(regimes))
  names(nodes) <- as.character(regimes)

  # start with all regimes being descending from the root; we will update them
  # as we encounter them in a preorder traversal.
  nodes[] <- 0L

  eFromRoot <- which(tree$edge[, 1] == N + 1L)
  rootRegime <- sort(tree$edge.regime[eFromRoot])[1L]
  nodes[as.character(rootRegime)] <- N + 1L

  for(ei in preorder) {
    if(as.character(tree$edge.regime[ei]) %in% regimes &&
       nodes[as.character(tree$edge.regime[ei])] == 0L) {
      i <- tree$edge[ei, 2L]
      j <- tree$edge[ei, 1L]

      nodes[as.character(tree$edge.regime[ei])] <- i
    }
  }

  nodes
}

#' Get the regimes of the branches leading to a set of nodes or tips
#' @param tree a phylo object with an edge.regime member denoting regimes.
#' @param nodes an integer vector denoting the nodes
#' @return a character or an integer vector denoting the regimes of the branches
#' leading to the nodes, according to tree$edge.regime.
#' @export
PCMTreeGetRegimesForNodes <- function(tree, nodes) {
  if(!inherits(tree, "phylo")) {
    stop("ERR:026h0:PCMBase:PCMTree.R:PCMTreeGetRegimesForNodes:: argument tree should be a phylo.")
  }

  tree$edge.regime[match(nodes, tree$edge[, 2])]
}

#' Get the tips belonging to a regime tree
#' @param tree a phylo object with an edge.regime member
#' @param regime a character or integer belonging to tree$edge.regime
#' @return an integer vector with the ids of the tips belonging to regime
#' @export
PCMTreeGetTipsInRegime <- function(tree, regime) {
  N <- PCMTreeNumTips(tree)
  tipEdges <- tree$edge[, 2] <= N & tree$edge.regime == regime
  tree$edge[tipEdges, 2]
}

#' Unique regimes on a tree in the order of occurrence from the root to the tips (preorder)
#'
#' @param tree a phylo object with an additional member edge.regime which should
#' be a character or an integer vector of length equal to the number of branches.
#' @param preorder an integer vector of row-indices in tree$edge matrix as returned
#' by PCMTreePreorder. This can be given for performance speed-up when several
#' operations needing preorder are executed on the tree. Default : \code{PCMTreePreorder(tree)}.
#' @return a character or an integer vector depending on tree$edge.regime.
#' @export
PCMTreeUniqueRegimes <- function(tree, preorder = PCMTreePreorder(tree)) {
  if(is.null(tree$edge.regime)) {
    stop("ERR:02610:PCMBase:PCMTree.R:PCMTreeUniqueRegimes:: tree$edge.regime is NULL,
         but should be a character or an integer vector denoting regime names.")
  }
  uniqueRegimesPos <- integer(PCMTreeNumUniqueRegimes(tree))

  uniqueRegimes <- sort(unique(tree$edge.regime))


  N <- PCMTreeNumTips(tree)
  eFromRoot <- which(tree$edge[, 1] == N+1)
  rootRegime <- sort(tree$edge.regime[eFromRoot])[1]

  names(uniqueRegimesPos) <- as.character(uniqueRegimes)

  uniqueRegimesPos[] <- 0L
  uniqueRegimesPos[as.character(rootRegime)] <- 1

  currentPos <- 2
  for(ei in preorder) {
    if(uniqueRegimesPos[as.character(tree$edge.regime[ei])] == 0L) {
      uniqueRegimesPos[as.character(tree$edge.regime[ei])] <- currentPos
      currentPos <- currentPos + 1
    }
  }

  uniqueRegimes[order(uniqueRegimesPos)]
}

#' Number of unique regimes on a tree
#' @param tree a phylo object
#' @return the number of different regimes encountered on the tree branches
#' @export
PCMTreeNumUniqueRegimes <- function(tree) {
  length(unique(tree$edge.regime))
}

#' Jumps in modeled traits associated with branches in a tree
#' @inheritParams PCMTreeNumTips
#' @return an integer vector of 0's and 1's with entries correspondin to the
#' denoting if a jump took place at the beginning of a branch.
#' @export
PCMTreeJumps <- function(tree) {
  if(!is.null(tree$edge.jump)) {
    if(length(tree$edge.jump) != nrow(tree$edge) ||
       !isTRUE(all(tree$edge.jump %in% c(0L, 1L)))) {
      stop("ERR:02681:PCMBase:PCMTree.R:PCMTreeJumps:: tree$edge.jump should
           be an integer vector of 0's and 1's with with length equal to the
           number of rows in tree$edge.")
    }
    as.integer(tree$edge.jump)
    } else {
      rep(0L, nrow(tree$edge))
    }
  }

#' Regimes associated with branches in a tree
#' @param tree a phylo object
#' @param model a PCM object
#' @param preorder an integer vector of row-indices in tree$edge matrix as returned
#' by PCMTreePreorder. This can be given for performance speed-up when several
#' operations needing preorder are executed on the tree. Default : \code{PCMTreePreorder(tree)}.
#' @return an integer vector with entries corresponding to the rows in tree$edge
#'   denoting the regime on each branch in the tree as an index in PCMRegimes(model).
PCMTreeMatchRegimesWithModel <- function(tree, model, preorder = PCMTreePreorder(tree)) {
  if(is.null(tree$edge.regime)) {
    PCMTreeSetDefaultRegime(tree, model)
  }
  #regimes <- match(tree$edge.regime, PCMRegimes(model, tree, preorder))
  regimes <- match(tree$edge.regime, PCMRegimes(model))
  if(any(is.na(regimes))) {
    stop(paste0("ERR:02671:PCMBase:PCMTree.R:PCMTreeMatchRegimesWithModel:: ",
                " Some of the regimes in tree$edge.regime not found in",
                "attr(model, 'regimes').\n",
                "unique regimes on the tree:", toString(PCMTreeUniqueRegimes(tree, preorder)), "\n",
                "attr(model, 'regimes')", toString(PCMRegimes(model))))
  }
  regimes
}

#' Pre-order tree traversal
#' @param tree a phylo object with possible singleton nodes (i.e. internal nodes with
#' one daughter node)
#' @return a vector of indices of edges in tree$edge in pre-order.
#' @export
PCMTreePreorder <- function(tree) {
  # number of tips
  N <- length(tree$tip.label)

  # total number of nodes in the tree is the number of edges + 1 for the root
  M <- dim(tree$edge)[1]+1

  ordFrom <- order(tree$edge[,1])

  # we need the ordered edges in order to easily traverse all edges starting from
  # a given node
  iFrom <- match(1:M, tree$edge[ordFrom, 1])

  # the result is a vector of edge indices in the breadth-first search order
  res <- vector(mode='integer', length=M-1)

  # node-indices at the current level (start from the root)
  cn <- N+1
  j <- 1
  while(length(cn)>0) {
    cnNext <- c()
    for(n in cn) {
      # if not a tip
      if(n > N) {
        # indices in ordFrom of edges starting from the current node
        if(n < M) {
          es <- iFrom[n]:(iFrom[n+1]-1)
        } else {
          es <- iFrom[n]:(M-1)
        }
        jNext <- j+length(es)
        res[j:(jNext-1)] <- ordFrom[es]
        j <- jNext
        cnNext <- c(cnNext, tree$edge[ordFrom[es], 2])
      }
    }
    cn <- cnNext
  }
  res
}

#' Post-order tree traversal
#' @param tree a phylo object with possible singleton nodes (i.e. internal nodes with
#' one daughter node)
#' @return a vector of indices of edges in tree$edge in post-order.
#' @export
PCMTreePostorder <- function(tree) {
  rev(PCMTreePreorder(tree))
}

#' A matrix (table) of ancestors/descendants for each node in a tree
#' @details This function has time and memory complexity O(M^2), where M is the
#'  number of nodes in the tree. It can take several minutes and gigabytes of
#'  memory on trees of more than 10000 tips.
#' @param tree a phylo object
#' @param preorder an integer vector returned by a previous call to \code{PCMTreePreorder(tree)}. Default \code{PCMTreePreorder(tree)}.
#' @return an integer square matrix of size M x M where M is the number of nodes
#' in the tree. Element j on row i is 0 if j is not an ancestor of i or a positive
#' integer equal to the position of j on the path from the root to i if j is an
#' ancestor of i.
#' @export
PCMTreeTableAncestors <- function(tree, preorder = PCMTreePreorder(tree)) {
  M <- PCMTreeNumNodes(tree)
  nodeLabels <- PCMTreeGetLabels(tree)

  tableAncestors <- matrix(0L, M, M)

  for(ei in preorder) {
    i <- tree$edge[ei, 2]
    j <- tree$edge[ei, 1]

    tableAncestors[i, ] <- tableAncestors[j, ]
    tableAncestors[i, j] <- max(tableAncestors[i,]) + 1
  }

  rownames(tableAncestors) <- colnames(tableAncestors) <- nodeLabels

  tableAncestors
}

#' Calculate the time from the root to each node of the tree
#' @param tree an object of class phylo
#' @param tipsOnly Logical indicating whether the returned results should be truncated only to the tips of the tree.
#' @return A vector of size the number of nodes in the tree (tips, root, internal) containing the time from the root to the corresponding node in the tree.
#' @export
PCMTreeNodeTimes <- function(tree, tipsOnly=FALSE) {
  preorder <- PCMTreePreorder(tree)
  es <- tree$edge[preorder, ]
  nEdges <- dim(es)[1]
  ts <- tree$edge.length[preorder]
  nodeTimes <- rep(0, PCMTreeNumNodes(tree))
  for(e in 1:nEdges)
    nodeTimes[es[e, 2]] <- nodeTimes[es[e, 1]]+ts[e]
  if(tipsOnly) {
    nodeTimes[1:PCMTreeNumTips(tree)]
  } else {
    nodeTimes
  }
}

#' A vector of the daughter nodes for a given parent node id in a tree
#' @param tree a phylo object.
#' @param parentId an integer denoting the id of the parent node
#' @return an integer vector of the direct descendants of parentId
#' @export
PCMTreeGetDaughters <- function(tree, parentId) {
  if(is.character(parentId)) {
    stop(paste0("ERR:026k1:PCMBase:PCMTree.R:PCMTreeGetParent::",
                " parentId should be integer but was character"))
  } else {
    tree$edge[tree$edge[, 1] == parentId, 2]
  }
}

#' The parent node id of a daughter node in a tree
#' @param tree a phylo object.
#' @param daughterId an integer denoting the id of the daughter node
#' @return an integer denoting the parent of daughterId
#' @export
PCMTreeGetParent <- function(tree, daughterId) {
  if(is.character(daughterId)) {
    stop(paste0("ERR:026j1:PCMBase:PCMTree.R:PCMTreeGetParent::",
                     " daughterId should be integer but was character"))
  } else {
    tree$edge[tree$edge[, 2] == daughterId, 1]
  }

}

#' The length of the branch leading to a node
#' @param tree a phylo object.
#' @param daughterId an integer denoting the id of a daughter node
#' @return a double denoting the length of the branch leading to daughterId
#' @export
PCMTreeGetBranchLength <- function(tree, daughterId) {
  tree$edge.length[tree$edge[, 2] == daughterId]
}

#' A list of the descendants for each node in a tree
#' @details This function has time and memory complexity O(M^2), where M is the
#'  number of nodes in the tree. It can take several minutes and gigabytes of
#'  memory on trees of more than 10000 tips.
#' @param tree a phylo object
#' @param tableAncestors an integer matrix resulting from a call to
#' PCMTreeTableAncestors(tree).
#' @return a list with unnamed elements in the order of nodes in the tree. Each
#' element is an integer vector containing the descendant nodes (in increasing
#'  order) of the node identified by its index-number in the list.
#' @export
PCMTreeListDescendants <- function(tree, tableAncestors = PCMTreeTableAncestors(tree)) {
  apply(tableAncestors, 2, function(descj) which(descj>0))
}

#' A list of the path to the root from each node in a tree
#' @details This function has time and memory complexity O(M^2), where M is the
#'  number of nodes in the tree. It can take several minutes and gigabytes of
#'  memory on trees of more than 10000 tips.
#' @param tree a phylo object
#' @param tableAncestors an integer matrix resulting from a call to
#' PCMTreeTableAncestors(tree).
#' @return a list with unnamed elements in the order of nodes in the tree. Each
#' element is an integer vector containing the ancestors nodes on the path from
#' the node (i) to the root of the tree in that order (the first element in the
#'  vector is the parent node of i and so on).
#' @export
PCMTreeListRootPaths <- function(tree, tableAncestors = PCMTreeTableAncestors(tree)) {
  apply(tableAncestors, 1, function(anci) {
    path <- which(anci>0)
    path[order(anci[path], decreasing = TRUE)]
  })
}

#' A list of all possible clade partitions of a tree with a number of splitting nodes
#'
#' @details Each subset of \code{nNodes} distinct internal or tip nodes
#' defines a partitioning of the branches of the tree into \code{nNodes+1} blocks.
#' This function generates partitions in which \code{nNode} of the blocks
#' are monophyletic complete groups (clades), while the \code{(nNodes+1)}'th block
#' is a subtree originating at the root with tips ending at the rooting nodes of
#' the \code{nNode} clades, eventually containing a clade of tips.
#' @param tree a phylo object
#' @param nNodes an integer giving the number of partitioning nodes. There would be
#' \code{nNodes+1} blocks in each partition (see details).
#' @param minCladeSize integer indicating the minimum number of tips allowed in a clade.
#' @param skipNodes integer vector indicating nodes that should not be used as
#' partition nodes. By default, this is an empty vector.
#' @param tableAncestors NULL (default) or an integer matrix returned by a previous call
#' to \code{PCMTreeTableAncestors(tree)}.
#' @param verbose a logical indicating if informative messages should be printed to
#' the console.
#'
#' @return a list of integer \code{nNodes}-vectors.
#' @export
PCMTreeListCladePartitions <- function(
  tree, nNodes, minCladeSize = 0, skipNodes = integer(0),
  tableAncestors = NULL, verbose = FALSE) {
  envir <- new.env()

  envir$M <- PCMTreeNumNodes(tree)
  envir$N <- PCMTreeNumTips(tree)
  envir$nNodes <- nNodes
  envir$minCladeSize <- minCladeSize

  if(!is.null(tableAncestors)) {
    envir$tableAncestors <- tableAncestors
  } else {
    envir$tableAncestors <- PCMTreeTableAncestors(tree)
  }
  envir$listDesc <- PCMTreeListDescendants(tree, envir$tableAncestors)
  envir$listRootPaths <- PCMTreeListRootPaths(tree, envir$tableAncestors)


  envir$nodesInUse <- rep(FALSE, envir$M)
  envir$nodesInUse[skipNodes] <- TRUE
  envir$nodesInUse[envir$N+1] <- TRUE
  envir$numNodesInUse <- sum(envir$nodesInUse)

  envir$listPartitions <- list()
  envir$nextPartition <- 1L
  envir$numTries <- 0L
  envir$numRemainingTips <- envir$N

  addToPartition <- function(partition, i, envir) {

    numTips <- sum(envir$listDesc[[i]] <= envir$N)
    if(envir$nodesInUse[i] ||
       numTips < envir$minCladeSize ||
       (envir$numRemainingTips - numTips) < envir$minCladeSize) {
      envir$numTries <- envir$numTries + 1L
      return(NULL)
    } else {
      partition <- c(partition, i)

      # ATTENTION! some of the elements of envir$listRootPaths[[i]] might
      # already be in use.
      newNodesInUse <- c(
        envir$listRootPaths[[i]][ !envir$nodesInUse[envir$listRootPaths[[i]]] ],
        i,
        envir$listDesc[[i]][ !envir$nodesInUse[envir$listDesc[[i]]] ])

      numNewNodesInUse <- length(newNodesInUse)

      envir$nodesInUse[newNodesInUse] <- TRUE
      envir$numNodesInUse <- envir$numNodesInUse + numNewNodesInUse
      envir$numRemainingTips <- envir$numRemainingTips - numTips

      if(length(partition) == envir$nNodes) {
        envir$numTries <- envir$numTries + 1L

        envir$listPartitions[[envir$nextPartition]] <- partition
        if(verbose && envir$nextPartition %% 1000 == 0) {
          cat("Generated ", envir$nextPartition, " partitions out of ",
              envir$numTries, " ...\n")
        }

        envir$nextPartition <- envir$nextPartition + 1L
      } else {
        if(envir$numNodesInUse < envir$M) {
          for(iNext in (i+1):envir$M) {
            addToPartition(partition, iNext, envir)
          }
        }
      }
      envir$nodesInUse[newNodesInUse] <- FALSE
      envir$numNodesInUse <- envir$numNodesInUse - numNewNodesInUse
      envir$numRemainingTips <- envir$numRemainingTips + numTips
    }
  }

  for(i in 1:envir$M) {
    if(verbose) {
      cat("Trying with first node in partition: ", i, "...\n")
    }
    addToPartition(c(), i, envir)
  }

  envir$listPartitions
}

#' A list of all possible (including recursive) partitions of a tree
#'
#' @param tree a phylo object with set labels for the internal nodes
#' @param minCladeSize integer indicating the minimum number of tips allowed in
#' one part.
#' @param tableAncestors NULL (default) or an integer matrix returned by a
#' previous call to \code{PCMTreeTableAncestors(tree)}.
#' @param verbose a logical indicating if informative messages should be printed to
#' the console.
#'
#' @return a list of integer vectors.
#' @export
PCMTreeListAllPartitions <- function(
  tree,
  minCladeSize,
  tableAncestors = NULL,
  verbose = FALSE) {

  if(is.null(tableAncestors)) {
    if(verbose) {
      cat("Creating tableAncestors...\n")
    }
    tableAncestors <- PCMTreeTableAncestors(tree)
  } else {
    colnames(tableAncestors) <- rownames(tableAncestors) <- PCMTreeGetLabels(tree)
  }

  res <- PCMTreeListAllPartitionsInternal(
    tree = tree,
    minCladeSize = minCladeSize,
    withoutNodesLabels = character(0),
    tableAncestors = tableAncestors,
    verbose = verbose,
    level = 0L
  )

  lapply(res, function(p) PCMTreeMatchLabels(tree, p))
}

PCMTreeListAllPartitionsInternal <- function(
  tree,
  minCladeSize,
  withoutNodesLabels,
  tableAncestors,
  verbose = FALSE,
  level = 0L) {

  if(verbose) {
    indent <- do.call(paste0, as.list(rep("  ", level)))
  }

  N <- PCMTreeNumTips(tree)
  M <- PCMTreeNumNodes(tree)
  labels <- PCMTreeGetLabels(tree)
  rootNodeLabel <- labels[N + 1L]

  if(verbose) {
    cat(indent,
        "PCMTreeListAllPartitionsInternal called on a tree with", N, "tips and", M,
        "nodes. Skipping ",
        length(intersect(PCMTreeGetLabels(tree), withoutNodesLabels)),
        "nodes\n")
  }

  partitionNodesLabels <- labels[
    unlist(PCMTreeListCladePartitions(
      tree = tree, nNodes = 1L, minCladeSize = minCladeSize,
      skipNodes = PCMTreeMatchLabels(
        tree, intersect(PCMTreeGetLabels(tree), withoutNodesLabels)),
      tableAncestors = tableAncestors, verbose = FALSE))]

  if(verbose) {
    cat(indent, "partitionNodesLabels:\n")
    print(partitionNodesLabels)
  }

  if(length(partitionNodesLabels) == 0L) {
    list(character(0))
  } else {
    iLabel <- partitionNodesLabels[1]

    # The set of all partitions of tree can be divided in two subsets:
    # 1. the subset containign all partitions without node i
    # 2. ths subset containing all partitions with node i

    if(verbose) {
      cat(indent, "1. find all partitions without node", iLabel, "\n")
    }
    partitionsWithouti <- PCMTreeListAllPartitionsInternal(
      tree = tree,
      minCladeSize = minCladeSize,
      withoutNodesLabels = c(withoutNodesLabels, iLabel),
      tableAncestors = tableAncestors,
      verbose = verbose,
      level = level + 1L)

    if(verbose) {
      print(partitionsWithouti)
    }

    if(verbose) {
      cat(indent, "2. find all partitions with node", iLabel, "\n")
    }

    # 2. list all partitions with node i:
    if(verbose) {
      cat(indent, "Splitting tree at node", iLabel, "\n")
    }
    spl <- PCMTreeSplitAtNode(
      tree = tree, node = iLabel, tableAncestors = tableAncestors)

    partitionsClade <- PCMTreeListAllPartitionsInternal(
      tree = spl$clade,
      minCladeSize = minCladeSize,
      withoutNodesLabels = withoutNodesLabels,
      tableAncestors = tableAncestors[
        PCMTreeGetLabels(spl$clade), PCMTreeGetLabels(spl$clade)],
      verbose = verbose,
      level = level + 1L)

    if(verbose) {
      cat(indent, "partitionsClade:\n")
      print(partitionsClade)
    }

    partitionsRest <- PCMTreeListAllPartitionsInternal(
      tree = spl$rest,
      minCladeSize = minCladeSize,
      withoutNodesLabels = withoutNodesLabels,
      tableAncestors = tableAncestors[
        PCMTreeGetLabels(spl$rest), PCMTreeGetLabels(spl$rest)],
      verbose = verbose,
      level = level + 1L)

    res <- partitionsWithouti

    for(k in seq_along(partitionsClade)) {
      for(l in seq_along(partitionsRest)) {
        res[[length(res) + 1L]] <-
          c(partitionsRest[[l]], iLabel, partitionsClade[[k]])
        if(verbose) {
          cat(indent, "Adding\n")
          print(c(partitionsRest[[l]], iLabel, partitionsClade[[k]]))
        }
      }
    }
    res
  }
}


#' Slit a tree at a given internal node into a clade rooted at this node and the remaining tree after dropping this clade
#' @param tree a phylo object
#' @param node an integer or character indicating a root, internal or tip node
#' @param tableAncestors an integer matrix returned by a previous call to PCMTreeTableAncestors(tree) or NULL.
#' @param X an optional k x N matrix with trait value vectors for each tip in tree.
#' @return A list containing two named phylo objects:
#' \itemize{
#' \item{clade }{The subtree (clade) starting at \code{node}.}
#' \item{Xclade }{The portion of X attributable to the tips in clade; NULL if X is NULL.}
#' \item{rest }{The tree resulting after dropping all tips in the clade.}
#' \item{Xrest }{The portion of X attributable to the tips in rest; NULL if X is NULL.}
#' }
#' @details In the current implementation, the edge.jump and edge.regime members
#' of the tree will be discarded and not present in the clade.
#'
#' TODO: preserveRegimes and preserve regimes
#'
#' @importFrom ape drop.tip
#' @export
PCMTreeSplitAtNode <- function(tree, node, tableAncestors = PCMTreeTableAncestors(tree), X=NULL) {
  if(!inherits(tree, "phylo")) {
    stop("ERR:02600:PCMBase:PCMTree.R:PCMTreeSplit:: tree must be a phylo object.")
  }
  N <- PCMTreeNumTips(tree)
  M <- PCMTreeNumNodes(tree)
  nodeLabels <- PCMTreeGetLabels(tree)
  if(is.character(node)) {
    nodeLab <- node
    node <- match(nodeLab, nodeLabels)
    if(is.na(node)) {
      stop(paste0("ERR:02601:PCMBase:PCMTree.R:PCMTreeSplit:: character node (", node, ") was not matched against the nodeLabels in tree."))
    }
  } else {
    nodeLab <- try(nodeLabels[node], silent = TRUE)
    if(class(nodeLab) == "try-error") {
      stop(paste0("ERR:02602:PCMBase:PCMTree.R:PCMTreeSplit:: non-character type node (", node, ") was not found: ", nodeLab, "."))
    }
  }

  # remove edge.regime and edge.jump from the tree - we currently do not update them.
  if( !is.null(tree$edge.regime) ) {
    tree <- tree[- which(names(tree) == "edge.regime")]
    class(tree) <- "phylo"
  }
  if(!is.null(tree$edge.jump)) {
    tree <- tree[- which(names(tree) == "edge.jump")]
    class(tree) <- "phylo"
  }

  if(node == N+1) {
    list(clade = tree,
         Xclade = X,
         rest = NULL,
         Xrest = NULL)
  } else {
    if(is.null(tableAncestors)) {
      tableAncestors <- PCMTreeTableAncestors(tree)
    } else {
      tableAncestors <- tableAncestors[nodeLabels, nodeLabels]
    }
    nodesClade <- which(tableAncestors[, node] > 0)
    tipsClade <- nodesClade[nodesClade <= N]

    if(length(tipsClade) == 0 && node <= N) {
      tipsClade <- node
    }

    if(!is.null(X)) {
      colnames(X) <- tree$tip.label <- as.character(1:N)
    }


    clade = drop.tip(tree, tip = setdiff(1:N, tipsClade), trim.internal = TRUE, collapse.singles = FALSE)

    # Because we used collapse.singles=FALSE, clade is still holding the entire lineage from the root to
    # node. We need to cut that part manually so that clade is rooted at node. We don't have to do this for rest,
    # because it should indeed be rooted at the tree-root.
    cladeNodeLabels <- PCMTreeGetLabels(clade)
    cladeEdgeNodeLabels <- cbind(cladeNodeLabels[clade$edge[, 1]], cladeNodeLabels[clade$edge[, 2]])
    cladeRootLab <- cladeNodeLabels[PCMTreeNumTips(clade) + 1]
    # TODO:  avoid this while loop and remove all root-nodes at once
    while(cladeRootLab != nodeLab) {
      ei <- which(cladeEdgeNodeLabels[, 1] == cladeRootLab)

      cladeRootLabNew <- cladeEdgeNodeLabels[ei, 2]
      clade$edge <- clade$edge[-ei, ]
      clade$edge[clade$edge > length(tipsClade)] <- clade$edge[clade$edge > length(tipsClade)] - 1
      clade$edge.length <- clade$edge.length[-ei]
      cladeEdgeNodeLabels <- cladeEdgeNodeLabels[-ei, ]
      clade$node.label <- clade$node.label[-match(cladeRootLab, clade$node.label)]
      clade$Nnode <- clade$Nnode - 1
      cladeRootLab <- cladeRootLabNew
    }

    rest = drop.tip(tree, tip = tipsClade, trim.internal = TRUE, collapse.singles = FALSE)
    if(!is.null(X)) {
      Xclade <- X[, clade$tip.label, drop = FALSE]
      Xrest <- X[, rest$tip.label, drop = FALSE]
    } else {
      Xclade <- NULL
      Xrest <- NULL
    }
    list(clade = clade,
         Xclade = Xclade,
         rest = rest,
         Xrest = Xrest)
  }
}


#' Extract a clade from phylogenetic tree
#' @param tree a phylo object
#' @param cladeRootNode a character string denoting the label or an integer denoting a node in the tree.
#' @param tableAncestors an integer matrix returned by a previous call to PCMTreeTableAncestors(tree) or NULL.
#' @param X an optional k x N matrix with trait value vectors for each tip in tree.
#' @param returnPhylo logical indicating if only the phylo object associated with the clade should be returned.
#' Defaults to \code{is.null(X)}
#' @return If returnPhylo is TRUE, a phylo object associated with the clade, otherise, a list with two named members :
#' \itemize{
#' \item{tree}{the phylo object associated with the clade}
#' \item{X}{the submatrix of X with columns corresponding to the tips in the clade}
#' }
#' @seealso PCMTreeSpliAtNode PCMTreeDropClade
#' @export
PCMTreeExtractClade <- function(tree, cladeRootNode, tableAncestors = NULL, X=NULL, returnPhylo=is.null(X)) {
  if(is.character(cladeRootNode)) {
    if(!is.character(tree$node.label)) {
      stop(paste0("ERR:02620:PCMBase:PCMTree.R:PCMTreeClade:", cladeRootNode, ": cladeRootNode is a character string but tree$node.label is missing or not a character vector."))
    } else {
      whichNode <- which(tree$node.label == cladeRootNode)
      whichTip <- which(tree$tip.label == cladeRootNode)
      if(length(whichNode) > 0) {
        cladeRootNodeNumber <- PCMTreeNumTips(tree) + whichNode
      } else if(length(whichTip) > 0) {
        cladeRootNodeNumber <- whichTip
      } else {
        cladeRootNodeNumber <- NA
      }

      if(is.na(cladeRootNodeNumber)) {
        stop(paste0("ERR:02621:PCMBase:PCMTree.R:PCMTreeClade:", cladeRootNode, ": cladeRootNode of character-type was not found in tree$node.label"))
      }
    }
  } else {
    cladeRootNodeNumber <- as.integer(cladeRootNode)
    if(cladeRootNodeNumber <= 0 || cladeRootNodeNumber > PCMTreeNumNodes(tree)) {
      stop(paste0("ERR:02622:PCMBase:PCMTree.R:PCMTreeClade:", cladeRootNode, ": cladeRootNode of integer should be between 1 and M=", PCMTreeNumNodes(tree), " (the number of nodes in the tree)."))
    }
  }

  spl <- PCMTreeSplitAtNode(tree, cladeRootNodeNumber, tableAncestors = tableAncestors, X = X)

  if(returnPhylo) {
    spl$clade
  } else {
    list(tree = spl$clade, X = spl$Xclade)
  }
}


#' Drop a clade from a phylogenetic tree
#' @param tree a phylo object
#' @param cladeRootNode a character string denoting the label or an integer denoting a node in the tree
#' @param tableAncestors an integer matrix returned by a previous call to PCMTreeTableAncestors(tree) or NULL.
#' @param X an optional k x N matrix with trait value vectors for each tip in tree.
#' @param returnPhylo logical indicating if only the phylo object associated with the tree after dropping the clade should be returned. Defaults to \code{is.null(X)}
#' @param errorOnMissing logical indicating if an error should be rased if cladeRootNode is not among the
#' nodes in tree. Default FALSE, meaning that if cladeRootNode is not a node in tree the tree (and X if
#' returnPhylo is FALSE) is/are returned unchanged.
#' @return If returnPhylo is TRUE, a phylo object associated with the remaining tree after dropping the clade, otherise, a list with two named members :
#' \itemize{
#' \item{tree}{the phylo object associated with the remaining tree after dropping the clade}
#' \item{X}{the submatrix of X with columns corresponding to the tips in the remaining tree}
#' }
#' @seealso PCMTreeSpliAtNode PCMTreeDropClade
#' @export
PCMTreeDropClade <- function(tree, cladeRootNode, tableAncestors = NULL, X=NULL, returnPhylo=is.null(X), errorOnMissing = FALSE) {
  if(is.character(cladeRootNode)) {
    if(!is.character(tree$node.label)) {
      stop(paste0("ERR:02630:PCMBase:PCMTree.R:PCMTreeClade:", cladeRootNode, ": cladeRootNode is a character string but tree$node.label is missing or not a character vector."))
    } else {
      whichNode <- which(tree$node.label == cladeRootNode)
      whichTip <- which(tree$tip.label == cladeRootNode)
      if(length(whichNode) > 0) {
        cladeRootNodeNumber <- PCMTreeNumTips(tree) + whichNode
      } else if(length(whichTip) > 0) {
        cladeRootNodeNumber <- whichTip
      } else {
        cladeRootNodeNumber <- NA
      }
    }
  } else {
    cladeRootNodeNumber <- as.integer(cladeRootNode)
  }

  res <- if(returnPhylo) {
    tree
  } else {
    list(tree = tree, X = X)
  }

  err <- NULL
  skipSplit <- FALSE
  if(is.na(cladeRootNodeNumber)) {
    skipSplit <- TRUE
    if(errorOnMissing) {
      err <- paste0("ERR:02631:PCMBase:PCMTree.R:PCMTreeClade:", cladeRootNode, ": cladeRootNode of character-type was not found in tree$node.label")
    }
  } else if(cladeRootNodeNumber <= 0 || cladeRootNodeNumber > PCMTreeNumNodes(tree)) {
    skipSplit <- TRUE
    if(errorOnMissing) {
      err <- paste0("ERR:02632:PCMBase:PCMTree.R:PCMTreeClade:", cladeRootNode, ": cladeRootNode of integer should be between 1 and M=", PCMTreeNumNodes(tree), " (the number of nodes in the tree).")
    }
  }

  if(!is.null(err)) {
    stop(err)
  }

  if(!skipSplit) {
    spl <- PCMTreeSplitAtNode(tree, cladeRootNodeNumber, tableAncestors = tableAncestors, X = X)

    if(returnPhylo) {
      res <- spl$rest
    } else {
      res <- list(tree = spl$rest, X = spl$Xrest)
    }
  }

  res
}


#' Perfrorm nested extractions or drops of clades from a tree
#' @param tree a phylo object with named tips and internal nodes
#' @param expr a character string representing an R expression of nested calls of functions
#' \code{E(x,node)} denoting extracting the clade rooted at node from the tree x, or \code{D(x,node)},
#' denoting dropping the clade rooted at node from the tree x. These calls can be nested, i.e. x can
#' be either the symbol x (corresponding to the original tree passed as argument) or a nested call to
#' d or e.
#' @return the resulting phylo object from evaluating expr on tree.
#' @export
PCMTreeEvalNestedEDxOnTree <- function(expr, tree) {

  tableAncestors <- PCMTreeTableAncestors(tree)
  env <- new.env()

  env$E <- function(x,node) {
    PCMTreeExtractClade(tree = x, cladeRootNode = as.character(node), tableAncestors = tableAncestors)
  }
  env$D <- function(x,node) {
    PCMTreeDropClade(tree = x, cladeRootNode = as.character(node), tableAncestors = tableAncestors)
  }
  env$x <- tree
  eval(parse(text=expr), envir = env)
}


#' A matrix with the begin and end time from the root for each edge in tree
#' @param tree a phylo
#' @export
PCMTreeEdgeTimes <- function(tree) {
  nodeTimes <- PCMTreeNodeTimes(tree)

  # begin and end-time of each edge relative to the root
  edgeTimes <- matrix(0.0, nrow(tree$edge), 2)
  edgeTimes[, 1] <- nodeTimes[tree$edge[, 1]]
  edgeTimes[, 2] <- nodeTimes[tree$edge[, 2]]
  edgeTimes
}

#' Find the crossing points of an epoch-time with each lineage of a tree
#' @param tree a phylo
#' @param epoch a positive numeric ing tip-ward distance from the root
#' @return a named list with an integer vector element "nodes" denoting the ending nodes for each
#' branch crossing epoch and numeric vector element "positions" denoting the root-ward offset
#' from each node in nodes.
#' @export
PCMTreeLocateEpochOnBranches <- function(tree, epoch) {
  edgeTimes <- PCMTreeEdgeTimes(tree)
  where <- (edgeTimes[, 1] < epoch) & (epoch <= edgeTimes[, 2])
  nodes <- tree$edge[where, 2]
  positions <- edgeTimes[where, 2] - epoch
  list(nodes = nodes, positions = positions)
}

#' Find the middle point of each branch longer than a threshold
#' @param tree a phylo
#' @param threshold a positive numeric; only branches longer than threshold will be returned; Default 0.
#' @return a named list with an integer vector element "nodes" denoting the ending nodes for each
#' branch crossing epoch and numeric vector element "positions" denoting the root-ward offset
#' from each node in nodes.
#' @export
PCMTreeLocateMidpointsOnBranches <- function(tree, threshold = 0) {
  where <- tree$edge.length > threshold
  nodes <- tree$edge[where, 2]
  positions <- tree$edge.length[where]/2
  list(nodes = nodes, positions = positions)
}

#' Insert singleton nodes on chosen edges
#'
#' @param tree a phylo object
#' @param nodes an integer vector denoting the terminating nodes of the edges on which
#' a singleton node is to be inserted
#' @param positions a positive numeric vector of the same length as nodes denoting
#' the root-ward distances from nodes at which the singleton nodes should be
#' inserted.
#' @param epoch a numeric indicating a distance from the root at which a singleton
#' node should be inserted in all lineages that are alive at that time.
#' @param minLength a numeric indicating the minimum allowed branch-length after
#' dividing a branch by insertion of a singleton nodes. No singleton node is inserted
#' if this would result in a branch shorter than `minLength`. Note that this condition
#' is checked only in `PCMTreeInsertSingletonsAtEpoch`.
#'
#' @importFrom ape bind.tree drop.tip
#' @return a modified version of tree with inserted singleton nodes at the specified locations
#' @seealso \code{\link{PCMTreeEdgeTimes}} \code{\link{PCMTreeLocateEpochOnBranches}} \code{\link{PCMTreeLocateMidpointsOnBranches}}
#' @export
PCMTreeInsertSingletons <- function(tree, nodes, positions) {
  M <- PCMTreeNumNodes(tree)
  N <- PCMTreeNumTips(tree)
  tip.label <- tree$tip.label
  node.label <- tree$node.label
  edge.regime <- tree$edge.regime
  edge.length <- tree$edge.length

  PCMTreeSetLabels(tree, paste0("x_", 1:PCMTreeNumNodes(tree)))

  if( !is.null(node.label) ) {
    names(node.label) <- paste0("x_", (N+1):M)
  }

  edgeNames <- paste0("x_", tree$edge[, 2])
  if( !is.null(edge.regime) ) {
    names(edge.regime) <- edgeNames
  }
  if( !is.null(edge.length) ) {
    names(edge.length) <- edgeNames
  }

  edgeTimes <- PCMTreeEdgeTimes(tree)
  rownames(edgeTimes) <- edgeNames

  # which edges should be processed
  names(nodes) <- names(positions) <- PCMTreeGetLabels(tree)[nodes]
  edgesToCut <- tree$edge[, 2] %in% nodes
  edgesToCutNames <- edgeNames[edgesToCut]

  edge.regime.new <- rep("", length(edgesToCutNames))
  names(edge.regime.new) <- paste0(
    "i_", round(edgeTimes[edgesToCutNames, 2] - positions[edgesToCutNames], 2),
    "_", edgesToCutNames)

  node.label.new <- names(edge.regime.new)
  names(node.label.new) <- names(edge.regime.new)

  tipTree <- list(edge = matrix(c(2, 1), nrow = 1, ncol = 2),
                  tip.label = "y_1",
                  node.label = "y_2",
                  edge.length = c(1.0),
                  Nnode = 1)
  class(tipTree) <- "phylo"

  for(edgeName in edgesToCutNames) {
    # find the number of the ending node for the edge
    node <- N + match(edgeName, tree$node.label)
    if(is.na(node)) {
      # node should be a tip
      node <- match(edgeName, tree$tip.label)
    }
    #cat("edgeName=", edgeName, "node=", node, "\n")

    # find the position (root-ward offset from node) where to cut

    # add, then drop the tipTree using ape's functions
    tree <- bind.tree(tree, tipTree, node, positions[edgeName])
    tree <- drop.tip(tree, "y_1", collapse.singles = FALSE)

    # inserted edge
    edgeNameNew <- paste0(
      "i_",
      round(edgeTimes[edgeName, 2] - positions[edgeName], 2), "_", edgeName)

    tree$node.label[is.na(tree$node.label)] <- edgeNameNew

    if(!is.null(edge.regime)) {
      edge.regime.new[edgeNameNew] <- edge.regime[edgeName]
    }
  }

  # The edge leading to the new internal node take the same regime as ITS
  # DAUGHTER edge
  if(!is.null(edge.regime)) {
    edge.regime.new <- c(edge.regime, edge.regime.new)
    nodeNamesNew <- PCMTreeGetLabels(tree)
    tree$edge.regime <- edge.regime.new[nodeNamesNew[tree$edge[, 2]]]
  }

  # restore original node labels
  if(!is.null(node.label)) {
    node.label.new <- c(node.label, node.label.new)
    node.label <- node.label.new[tree$node.label]
    tree$node.label <- node.label
  }

  # restore original tip.label
  if(!is.null(tip.label)) {
    tree$tip.label <- tip.label
  }

  tree
}

#' @describeIn PCMTreeInsertSingletons
#'
#' @export
PCMTreeInsertSingletonsAtEpoch <- function(tree, epoch, minLength = 0.1) {
  nodeTimes <- PCMTreeNodeTimes(tree)

  points <- PCMTreeLocateEpochOnBranches(tree, epoch)

  idxPoints <- sapply(seq_along(points$nodes), function(i) {
    node <- points$nodes[i]
    nodeTime <- nodeTimes[points$nodes[i]]
    parentNode <- tree$edge[tree$edge[, 2] == points$nodes[i], 1]
    parentTime <- nodeTimes[parentNode]
    pos <- points$positions[i]

    (points$positions[i] > minLength) && (nodeTime - pos - parentTime) > minLength
    })

  if(length(idxPoints) > 0) {
    points$nodes <- points$nodes[idxPoints]
    points$positions <- points$positions[idxPoints]
    if(length(points$nodes) > 0) {
      PCMTreeInsertSingletons(tree, points$nodes, points$positions)
    } else {
      tree
    }
  } else {
    tree
  }

}

#' Find the nearest node to a given time from the root (epoch) on each lineage crossing this epoch
#' @param tree a phylo
#' @param epoch a positive numeric
#' @return an integer vector
#' @export
PCMTreeNearestNodesToEpoch <- function(tree, epoch) {
  nodeTimes <- PCMTreeNodeTimes(tree)
  points <- PCMTreeLocateEpochOnBranches(tree, epoch)
  sapply(seq_along(points$nodes), function(i) {
    node <- points$nodes[i]
    parentNode <- tree$edge[tree$edge[, 2] == points$nodes[i], 1]

    if(points$positions[i] > (nodeTimes[points$nodes[i]] - points$positions[i] - nodeTimes[parentNode])) {
      parentNode
    } else {
      node
    }
  })
}

#' A data.table of the tips with their assigned regime
#' @param tree a phylo object with node-labels and regimes
#' @importFrom data.table data.table
#' @export
PCMTreeDtNodeRegimes <- function(tree) {
  N <- PCMTreeNumTips(tree)
  nodeTimes <- PCMTreeNodeTimes(tree)
  nodeLabels <- PCMTreeGetLabels(tree)

  data.table(
    startNode = tree$edge[, 1], endNode = tree$edge[, 2],
    startNodeLab = nodeLabels[tree$edge[, 1]], endNodeLab = nodeLabels[tree$edge[, 2]],
    startTime = nodeTimes[tree$edge[, 1]], endTime = nodeTimes[tree$edge[, 2]],
    regime = tree$edge.regime)
}

#' Prune the tree leaving one tip for each regime
#' @param tree a phylo with set node-labels and regimes
#' @return a pruned version of tree
#'
#' @seealso PCMTreeSetLabels PCMTreeSetRegimes
#'
#' @importFrom data.table data.table
#' @importFrom ape drop.tip
#'
#' @export
PCMTreeExtractBackboneRegimes <- function(tree) {

  # Needed to pass the check.
  regime <- endNode <- endTime <- NULL

  N <- PCMTreeNumTips(tree)
  nodeTimes <- PCMTreeNodeTimes(tree)
  nodeLabels <- PCMTreeGetLabels(tree)

  regimeNodesTree <- PCMTreeGetStartingNodesRegimes(tree)

  dtBranchRegimes <- data.table(
    startNode = tree$edge[, 1], endNode = tree$edge[, 2],
    startTime = nodeTimes[tree$edge[, 1]], endTime = nodeTimes[tree$edge[, 2]],
    regime = tree$edge.regime)

  dtTipRegimes <- dtBranchRegimes[endNode <= N]
  dtTipRegimes <- dtTipRegimes[, list(endNode=endNode[which.max(endTime)]), keyby = regime]

  tipsToDrop <- setdiff(1:N, dtTipRegimes[, endNode])

  tree2 <- drop.tip(tree, tipsToDrop, collapse.singles = FALSE)

  PCMTreeSetRegimes(tree2, nodeLabels[regimeNodesTree][-1], names(regimeNodesTree))
  tree2
}

#' A character representation of a phylo object.
#'
#' @param tree a phylo object.
#' @param includeLengths logical. Default: FALSE.
#' @param includeStartingNodesRegimes logical. Default: FALSE.
#' @return a character string.
#' @export
PCMTreeToString <- function(tree, includeLengths = FALSE, includeStartingNodesRegimes = FALSE) {

  orderEdge <- order(tree$edge[, 2])
  nodeLabs <- PCMTreeGetLabels(tree)

  edgeLabelsOrdered <- cbind(nodeLabs[tree$edge[orderEdge, 1]],
                             nodeLabs[tree$edge[orderEdge, 2]])
  attributes(edgeLabelsOrdered) <- NULL
  unname(edgeLabelsOrdered)

  if(includeLengths) {
    edgeLengthsOrdered <- tree$edge.length[orderEdge]
    attributes(edgeLengthsOrdered) <- NULL
    unname(edgeLengthsOrdered)
  } else {
    edgeLengthsOrdered <- ""
  }
  if(includeStartingNodesRegimes) {
    startingNodesRegimes <- PCMTreeGetStartingNodesRegimes(tree)
    startingNodesRegimesLabels <- nodeLabs[startingNodesRegimes]
    attributes(startingNodesRegimesLabels) <- NULL
    unname(startingNodesRegimesLabels)
  } else {
    startingNodesRegimesLabels <- ""
  }

  paste0(toString(edgeLabelsOrdered), "; ",
                         toString(edgeLengthsOrdered), "; ",
                         toString(startingNodesRegimesLabels))
}


#' Which tips in a tree belong to the same regime?
#' @param tree a phylo object with an existing edge.regime
#' @param upperTriangle logical indicating if all duplicated entries and diagonal
#' entries sould be set to NA (by default TRUE).
#' @param returnVector logical indicating if a vector instead of a matrix should be
#' returned (corresponding to calling as.vector on the resulting matrix and removing
#' NAs). Default: TRUE
#' @return a N x N logical matrix with TRUE on the diagonal and for each couple of
#' tips that belong to the same dagonal. if returnVector is TRUE (default) only a
#' vector of the non-NA entries will be returned.
#' @export
PCMTreeMatrixTipsInSameRegime <- function(tree, upperTriangle = TRUE, returnVector = TRUE) {
  nMax <- PCMTreeNumTips(tree)
  mat <- matrix(FALSE, nMax, nMax)
  diag(mat) <- TRUE

  for(i in 1:nMax) {
    for(j in 1:nMax) {
      ei <- which(tree$edge[, 2] == i)
      ej <- which(tree$edge[, 2] == j)
      mat[i,j] <- mat[j,i] <- (tree$edge.regime[ei] == tree$edge.regime[ej])
    }
  }
  if(upperTriangle) {
    mat[lower.tri(mat, diag = TRUE)] <- NA
  }
  if(returnVector) {
    mat <- as.vector(mat)
    mat <- mat[!is.na(mat)]
  }
  mat
}

#' Which nodes in a tree belong to the same regime?
#' @param tree a phylo object with an existing edge.regime
#' @param upperTriangle logical indicating if all duplicated entries and diagonal
#' entries sould be set to NA (by default TRUE).
#' @param returnVector logical indicating if a vector instead of a matrix should be
#' returned (corresponding to calling as.vector on the resulting matrix and removing
#' NAs). Default: TRUE
#' @return a M x M logical matrix with TRUE on the diagonal and for each couple of
#' tips that belong to the same dagonal. if returnVector is TRUE (default) only a
#' vector of the non-NA entries will be returned.
#' @export
PCMTreeMatrixNodesInSameRegime <- function(tree, upperTriangle = TRUE, returnVector = TRUE) {
  N <- PCMTreeNumTips(tree)
  nMax <- PCMTreeNumNodes(tree)
  mat <- matrix(FALSE, nMax, nMax)
  diag(mat) <- TRUE

  for(i in 1:nMax) {
    for(j in 1:nMax) {
      if( (i == N + 1) || (j == N + 1) ) {
        next
      } else {
        ei <- which(tree$edge[, 2] == i)
        ej <- which(tree$edge[, 2] == j)
        mat[i,j] <- mat[j,i] <- (tree$edge.regime[ei] == tree$edge.regime[ej])
      }
    }
  }
  mat[N+1,] <- mat[,N+1] <- NA

  if(upperTriangle) {
    mat[lower.tri(mat, diag = TRUE)] <- NA
  } else {

  }
  if(returnVector) {
    mat <- as.vector(mat)
    mat <- mat[!is.na(mat)]
  }
  mat
}


#' Plot a tree with regimes
#' @param tree a phylo with set labels and regimes
#' @param palette a named vector of colors corresponding to the regimes in tree
#' @param ... Arguments passed to ggtree, e.g. layout = 'fan', open.angle = 8, size=.25.
#' @note Currently, the ggtree package is not on CRAN and therefore it is not explicitly
#' imported by PCMBase.
#' @importFrom data.table data.table
#' @importFrom grDevices hcl
#' @importFrom ggplot2 aes scale_color_manual
#' @export
PCMTreePlot <- function(tree, palette = PCMColorPalette(PCMTreeNumUniqueRegimes(tree), PCMTreeUniqueRegimes(tree)), ...) {
  # Needed to pass the check
  regime <- NULL

  if(requireNamespace("ggtree")) {
    N <- PCMTreeNumTips(tree)
    R <- PCMTreeNumUniqueRegimes(tree)

    data <- rbind(data.table(node = tree$edge[, 2],
                             regime = as.factor(tree$edge.regime)),
                  data.table(node = N+1, regime = as.factor(tree$edge.regime[1])))

    plotTree <- ggtree::`%<+%`(ggtree::ggtree(tree, ...), data)

    plotTree + aes(color = regime) +
      scale_color_manual(name = "regime", values = palette)
  } else {
    stop("ERR:026i0:PCMBase:PCMTree.R:PCMTreePlot:: Calling PCMTreePlot needs ggtree package to be installed from Bioconductor. Check the instructions at https://bioconductor.org/packages/release/bioc/html/ggtree.html. Ggtree was not on CRAN at the time of releasing PCMBase and is not declared as dependency in the PCMBase-description.")
  }
}
