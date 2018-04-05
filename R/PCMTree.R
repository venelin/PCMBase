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
  nrow(tree$edge) + 1
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
#' @param nodes an integer vector denoting tip or internal nodes in tree - the
#' regimes change at the start of the branches leading to these nodes.
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
    if(length(regimes) != length(nodes) + 1 ||
       length(unique(regimes)) < length(regimes)) {
      stop("ERR:026d1:PCMBase:PCMTree.R:PCMTreeSetRegimes:: regimes should be a character or integer vector of length equal to length(nodes) + 1 and not have duplicated elements.")
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

#' Get the starting branch' nodes for each regime on a tree
#'
#' @details We call a starting branch the first branch from the root to the tips
#' with a given regime. A starting node is the node at which a starting branch
#' ends.
#' @param tree a phylo object with an edge.regime member denoting regimes. The
#' function assumes that each regime covers a linked set of branches on the tree.
#' @return an integer with elements equal to the starting nodes for each regime in
#' \code{regimes}.
#' @seealso \code{\link{PCMTreeSetRegimes}}
#' @export
PCMTreeGetStartingNodesRegimes <- function(tree) {
  if(!inherits(tree, "phylo")) {
    stop("ERR:026g0:PCMBase:PCMTree.R:PCMTreeGetStartNodesRegimes:: argument tree should be a phylo.")
  }

  regimes <-  PCMTreeUniqueRegimes(tree)

  N <- PCMTreeNumTips(tree)
  nodes <- integer(length(regimes))
  names(nodes) <- as.character(regimes)

  # start with all regimes being descending from the root; we will update them
  # as we encounter them in a preorder traversal.
  nodes[] <- 0

  eFromRoot <- which(tree$edge[, 1] == N + 1)
  rootRegime <- sort(tree$edge.regime[eFromRoot])[1]
  nodes[as.character(rootRegime)] <- N + 1

  preorder <- PCMTreePreorder(tree)

  for(ei in preorder) {
    if(as.character(tree$edge.regime[ei]) %in% regimes &&
       nodes[as.character(tree$edge.regime[ei])] == 0) {
      i <- tree$edge[ei, 2]
      j <- tree$edge[ei, 1]

      nodes[as.character(tree$edge.regime[ei])] <- i
    }
  }

  nodes
  }

#' Unique regimes on a tree in the order of occurrence from the root to the tips (preorder)
#'
#' @param tree a phylo object with an additional member edge.regime which should
#' be a character or an integer vector of length equal to the number of branches.
#'
#' @return a character or an integer vector depending on tree$edge.regime.
#' @export
PCMTreeUniqueRegimes <- function(tree) {
  if(is.null(tree$edge.regime)) {
    stop("ERR:02610:PCMBase:PCMTree.R:PCMRegimesUniqueTree:: tree$edge.regime is NULL,
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

  preorder <- PCMTreePreorder(tree)
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
#' @return an integer vector with entries corresponding to the rows in tree$edge
#'   denoting the regime on each branch in the tree as an index in PCMRegimes(model).
#' @export
PCMTreeMatchRegimesWithModel <- function(tree, model) {
  regimes <- match(tree$edge.regime, PCMRegimes(model, tree))
  if(any(is.na(regimes))) {
    stop(paste0("ERR:02671:PCMBase:PCMTree.R:PCMTreeMatchRegimesWithModel:: ",
                " Some of the regimes in tree$edge.regime not found in",
                "attr(model, 'regimes').\n",
                "unique regimes on the tree:", toString(PCMTreeUniqueRegimes(tree)), "\n",
                "attr(model, 'regimes')", toString(PCMRegimes(model))))
  }
  regimes
}

#' Pre-order tree traversal
#' @param tree a phylo object with possible singleton nodes (i.e. internal nodes with
#' one daughter node)
#' @return a vector of indices of edges in tree$edge in pre-order.
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


#' Extract information for pruning a tree used as cache during likelihood
#' calculation
#'
#' @param tree a phylo object
#'
#' @details The function is very slow on strongly unbalanced trees due to the
#' slow vector append() operation in R.
#'
#' @return a list of objects used by the R implementation of PCMLik().
#' @export
PCMTreePruningOrder <- function(tree) {
  # number of tips
  N <- length(tree$tip.label)
  # number of all nodes
  M <- nrow(tree$edge)+1

  # order the edge-indices in increasing index of ending node
  endingAt <- order(rbind(tree$edge, c(0, N+1))[, 2])

  edge <- tree$edge

  # count the number of non-visited children for each internal node
  nvc <- rep(0, M)

  # indices of parent node for edges that haven't still been gone through
  # initially, this are all edges
  ee1 <- edge[, 1]

  while(length(ee1)) {
    # For every element of (N+1):M its index in ee1 or NA
    matchp <- match((N+1):M, ee1)
    matchp <- matchp[!is.na(matchp)]
    # add one unvisited chlidren to each parent node's nvc
    nvc[ee1[matchp]] <- nvc[ee1[matchp]] + 1
    # remove the edges we've just gone through
    ee1 <- ee1[-matchp]
  }

  # start from the edges leading to tips
  nodesVector <- c()
  nodesIndex <- c(0)

  unVector <- c()
  # pointers to last unique indices (un) in unVector
  unIndex <- c(0)

  # internal or tip- nodes to which we are currently pointing, i.e. we are at
  # their parent-nodes and we are about to process the brances leading to them.
  nodes <- 1:N

  while(nodes[1] != N+1) {
    nodesIndex <- c(nodesIndex, nodesIndex[length(nodesIndex)]+length(nodes))
    nodesVector <- c(nodesVector, nodes)

    # indices of edges that end at nodes
    es <- endingAt[nodes]
    nodes <- c()

    while(length(es)>0) {
      # unique index of every edge ending at some of the nodes
      un <- match(unique(edge[es, 1]), edge[es, 1])
      # add these indices to unVector
      unVector <- c(unVector, un)
      # index of the last element of the current un in unVector
      unIndex <- c(unIndex, unIndex[length(unIndex)]+length(un))
      nvc[edge[es[un], 1]] <- nvc[edge[es[un], 1]] - 1
      nodes <- c(nodes, edge[es[un][nvc[edge[es[un], 1]] == 0], 1])
      es <- es[-un]
    }
  }
  list(# all raws from edge, times t and regimes must be accessed using indices
    # from the edingAt vector.
    endingAt=endingAt,

    nodesVector=nodesVector,
    nodesIndex=nodesIndex,
    nLevels=length(nodesIndex)-1,
    unVector=unVector,
    unIndex=unIndex)
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
  tableAncestors <- matrix(0, M, M)

  for(ei in preorder) {
    i <- tree$edge[ei, 2]
    j <- tree$edge[ei, 1]

    tableAncestors[i, ] <- tableAncestors[j, ]
    tableAncestors[i, j] <- max(tableAncestors[i,]) + 1
  }
  tableAncestors
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

#' A list of all possible partitions of a tree with a number of splitting nodes
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
#' @return a list of integer \code{nNodes}-vectors.
#' @export
PCMTreeListCladePartitions <- function(tree, nNodes, minCladeSize = 0, verbose = FALSE) {
  envir <- new.env()

  envir$M <- PCMTreeNumNodes(tree)
  envir$N <- PCMTreeNumTips(tree)
  envir$nNodes <- nNodes
  envir$minCladeSize <- minCladeSize

  envir$tableAncestors <- PCMTreeTableAncestors(tree)
  envir$listDesc <- PCMTreeListDescendants(tree, envir$tableAncestors)
  envir$listRootPaths <- PCMTreeListRootPaths(tree, envir$tableAncestors)


  envir$nodesInUse <- rep(FALSE, envir$M)
  envir$nodesInUse[envir$N+1] <- TRUE
  envir$numNodesInUse <- 1

  envir$listPartitions <- list()
  envir$nextPartition <- 1
  envir$numTries <- 0

  addToPartition <- function(partition, i, envir) {

    if(envir$nodesInUse[i] ||
       sum(envir$listDesc[[i]] <= envir$N) < envir$minCladeSize) {
      envir$numTries <- envir$numTries + 1
      return(NULL)
    } else {
      partition <- c(partition, i)

      # ATTENTION! some of the elements of envir$listRootPaths[[i]] might
      # already be in use.
      newNodesInUse <- c(
        envir$listRootPaths[[i]][ !envir$nodesInUse[envir$listRootPaths[[i]]] ],
        i,
        envir$listDesc[[i]])

      numNewNodesInUse <- length(newNodesInUse)

      envir$nodesInUse[newNodesInUse] <- TRUE
      envir$numNodesInUse <- envir$numNodesInUse + numNewNodesInUse

      if(length(partition) == envir$nNodes) {
        envir$numTries <- envir$numTries + 1

        envir$listPartitions[[envir$nextPartition]] <- partition
        if(verbose && envir$nextPartition %% 1000 == 0) {
          cat("Generated ", envir$nextPartition, " partitions out of ",
              envir$numTries, " ...\n")
        }

        envir$nextPartition <- envir$nextPartition + 1
      } else {
        if(envir$numNodesInUse < envir$M) {
          for(iNext in (i+1):envir$M) {
            addToPartition(partition, iNext, envir)
          }
        }
      }
      envir$nodesInUse[newNodesInUse] <- FALSE
      envir$numNodesInUse <- envir$numNodesInUse - numNewNodesInUse
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

#' Slit a tree at a given internal node into a clade rooted at this node and the remaining tree after dropping this clade
#' @param tree a phylo object
#' @param node an integer indicating an internal node
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
#' @importFrom ape drop.tip
#' @export
PCMTreeSplitAtNode <- function(tree, node, tableAncestors = PCMTreeTableAncestors(tree), X=NULL) {
  if(!inherits(tree, "phylo")) {
    stop("ERR:02600:PCMBase:PCMTree.R:PCMTreeSplit:: tree must be a phylo object.")
  }
  N <- PCMTreeNumTips(tree)

  if(!is.null(tree$edge.regime)) {
    tree <- tree[- which(names(tree) == "edge.regime")]
    class(tree) <- "phylo"
  }
  if(!is.null(tree$edge.jump)) {
    tree <- tree[- which(names(tree) == "edge.jump")]
    class(tree) <- "phylo"
  }

  if(node <= N) {
    stop("ERR:02601:PCMBase:PCMTree.R:PCMTreeSplit:: cannot split at a tip; node should be an internal node.")
  } else if(node == N+1) {
    list(clade = tree,
         Xclade = X,
         rest = NULL,
         Xrest = NULL)
  } else {
    if(is.null(tableAncestors)) {
      tableAncestors <- PCMTreeTableAncestors(tree)
    }
    nodesClade <- which(tableAncestors[, node] > 0)
    tipsClade <- nodesClade[nodesClade <= N]

    if(!is.null(X)) {
      colnames(X) <- tree$tip.label <- as.character(1:N)
    }
    clade = drop.tip(tree, tip = setdiff(1:N, tipsClade), trim.internal = TRUE, collapse.singles = TRUE)
    rest = drop.tip(tree, tip = tipsClade, trim.internal = TRUE, collapse.singles = TRUE)
    if(!is.null(X)) {
      Xclade <- X[, clade$tip.label]
      Xrest <- X[, rest$tip.label]
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
#' @param cladeRootNode a character string denoting the label or an integer denoting an internal node in tree
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
      cladeRootNodeNumber <- PCMTreeNumTips(tree) + which(tree$node.label == cladeRootNode)
      if(is.na(cladeRootNodeNumber)) {
        stop(paste0("ERR:02621:PCMBase:PCMTree.R:PCMTreeClade:", cladeRootNode, ": cladeRootNode of character-type was not found in tree$node.label"))
      }
    }
  } else {
    cladeRootNodeNumber <- as.integer(cladeRootNode)
    if(cladeRootNodeNumber <= PCMTreeNumTips(tree) || cladeRootNodeNumber > PCMTreeNumNodes(tree)) {
      stop(paste0("ERR:02621:PCMBase:PCMTree.R:PCMTreeClade:", cladeRootNode, ": cladeRootNode of integer type was either a tip (<=N) or bigger than M (the number of nodes in the tree)."))
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
#' @param cladeRootNode a character string denoting the label or an integer denoting an internal node in tree
#' @param tableAncestors an integer matrix returned by a previous call to PCMTreeTableAncestors(tree) or NULL.
#' @param X an optional k x N matrix with trait value vectors for each tip in tree.
#' @param returnPhylo logical indicating if only the phylo object associated with the tree after dropping the clade should be returned. Defaults to \code{is.null(X)}
#' @return If returnPhylo is TRUE, a phylo object associated with the remaining tree after dropping the clade, otherise, a list with two named members :
#' \itemize{
#' \item{tree}{the phylo object associated with the remaining tree after dropping the clade}
#' \item{X}{the submatrix of X with columns corresponding to the tips in the remaining tree}
#' }
#' @seealso PCMTreeSpliAtNode PCMTreeDropClade
#' @export
PCMTreeDropClade <- function(tree, cladeRootNode, tableAncestors = NULL, X=NULL, returnPhylo=is.null(X)) {
  if(is.character(cladeRootNode)) {
    if(!is.character(tree$node.label)) {
      stop(paste0("ERR:02620:PCMBase:PCMTree.R:PCMTreeClade:", cladeRootNode, ": cladeRootNode is a character string but tree$node.label is missing or not a character vector."))
    } else {
      cladeRootNodeNumber <- PCMTreeNumTips(tree) + which(tree$node.label == cladeRootNode)
      if(is.na(cladeRootNodeNumber)) {
        stop(paste0("ERR:02621:PCMBase:PCMTree.R:PCMTreeClade:", cladeRootNode, ": cladeRootNode of character-type was not found in tree$node.label"))
      }
    }
  } else {
    cladeRootNodeNumber <- as.integer(cladeRootNode)
    if(cladeRootNodeNumber <= PCMTreeNumTips(tree) || cladeRootNodeNumber > PCMTreeNumNodes(tree)) {
      stop(paste0("ERR:02621:PCMBase:PCMTree.R:PCMTreeClade:", cladeRootNode, ": cladeRootNode of integer type was either a tip (<=N) or bigger than M (the number of nodes in the tree)."))
    }
  }

  spl <- PCMTreeSplitAtNode(tree, cladeRootNodeNumber, tableAncestors = tableAncestors, X = X)

  if(returnPhylo) {
    spl$rest
  } else {
    list(tree = spl$rest, X = spl$Xrest)
  }
}

#' Perfrorm nested extractions (E) or drops (D) of clades from a tree
#' @param expr an R expression of nested calls of functions
#' \code{E(x,node)} denoting extracting the clade rooted at node from the tree x, or \code{D(x,node)},
#' denoting dropping the clade rooted at node from the tree x. These calls can be nested, i.e. x can
#' be a symbol or r expression evaluating to a  phylo object in
#' the global or calling environment (corresponding to the original tree passed
#'  as argument) or a nested call to D or E.
#' @return the resulting phylo object from evaluating expr in the calling environment
#' @export
PCMTreeEvalNestedEDs <- function(expr) {
  E <- function(x,node) {
    PCMTreeExtractClade(x, as.character(node))
  }
  D <- function(x,node) {
    PCMTreeDropClade(x, as.character(node))
  }
  eval(substitute(quote(expr)), parent.frame())
}
