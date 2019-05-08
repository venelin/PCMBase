# Copyright 2016-2019 Venelin Mitov
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


#' Create a PCMTree object from a phylo object
#'
#' @description PCMTree is class that inherits from the class 'phylo' in the
#' R-package 'ape'. Thus, all the functions working on a phylo object would work
#' in the same way if they recieve as argument an object of class 'PCMTree'. A
#' PCMTree object has the following members in addition to the regular members
#' ('tip.label', 'node.label', 'edge', 'edge.length') found in a regular phylo
#' object:
#' \describe{
#' \item{edge.part }{a character vector having as many elements as there are
#' branches in the tree (corresponding to the rows in `tree$edge`). Each element
#' denotes the name of the part to which the corresponding branch belongs. A
#' part in the tree represents a connected subset of its nodes and the branches
#' leading to these nodes. A partition of the tree represents the splitting of
#' the tree into a number of parts. Visually, a partition can be represented as
#' a coloring of the tree, in which no color is assigned to more than one part.
#' In other words, if two branches in the tree are connected by the same color,
#' they either share a node, or all the branches on the path in the tree
#' connecting these two branches have the same color. Formally, we define a
#' partition of the tree as any set of nodes in the tree that includes the root.
#' Each node in this set defines a part as the set of its descendant nodes that
#' can be reached without traversing another partition node. We name
#' each part by the label of its most ancestral node, that is, the node in it,
#' which is closest to the root fo the tree. The value of edge.part for an edge
#' in the tree is the name of the part that contsin the node ot which the edge
#' is pointing.}
#' \item{part.regime }{a named vector of size the number of parts in the tree.
#' The names correspond to part-names whereas the values denote the ids or
#' character names of regimes in a PCM object.}
#' }
#' The constructor PCMTree() returns an object of call
#' @param tree a phylo object. If this is already a PCMTree object, a copy of
#' this object will be returned.
#'
#' @return an object of class PCMTree. This is a copy of the passed phylo object
#' which is guaranteed to have node.label, edge.part and a part.regime
#' entries set.
#'
#' @examples
#' tree <- ape::rtree(8)
#'
#' # the following four are NULLs
#' tree$node.label
#' tree$edge.part
#' tree$part.regime
#' tree$edge.regime
#'
#' # In previous version regimes were assigned directly to the edges via
#' # tree$edge.regime. This is supported but not recommended anymore:
#'
#' tree$edge.regime <- sample(
#'   letters[1:3], size = PCMTreeNumNodes(tree) - 1, replace = TRUE)
#'
#' tree.a <- PCMTree(tree)
#' PCMTreeGetLabels(tree.a)
#' tree.a$node.label
#' tree.a$edge.part
#' tree.a$part.regime
#'
#' # this is set to NULL - starting from PCMBase 1.2.9 all of the information
#' # for the regimes is stored in tree$edge.part and tree$part.regime.
#' tree.a$edge.regime
#'
#' PCMTreeGetPartition(tree.a)
#' PCMTreeGetPartNames(tree.a)
#' PCMTreeGetPartRegimes(tree.a)
#'
#' # let's see how the tree looks like
#' \donttest{
#' PCMTreePlot(tree.a) + ggtree::geom_nodelab() + ggtree::geom_tiplab()
#' }
#'
#' # This is the recommended way to set a partition on the tree
#' PCMTreeSetPartition(tree.a, c(10, 12))
#'
#' PCMTreeGetPartition(tree.a)
#' PCMTreeGetPartNames(tree.a)
#' PCMTreeGetPartRegimes(tree.a)
#'
#'
#' \donttest{
#' PCMTreePlot(tree.a) + ggtree::geom_nodelab() + ggtree::geom_tiplab()
#' }
#' PCMTreeGetPartsForNodes(tree.a, c(11, 15, 12))
#' PCMTreeGetPartsForNodes(tree.a, c("11", "15", "12"))
#'
#' PCMTreeSetPartRegimes(tree.a, c(`9` = 'a', `10` = 'b', `12` = 'c'))
#'
#' PCMTreeGetPartition(tree.a)
#' PCMTreeGetPartNames(tree.a)
#' PCMTreeGetPartRegimes(tree.a)
#'
#' \donttest{
#' PCMTreePlot(tree.a) + ggtree::geom_nodelab() + ggtree::geom_tiplab()
#' }
#' @export
PCMTree <- function(tree) {
  if(is.PCMTree(tree)) {
    tree
  } else {
    if(!inherits(tree, "phylo")) {
      stop("PCMTree:: tree should be a phylo or a PCMTree object.")
    }

    class(tree) <- c("PCMTree", class(tree))

    if(is.null(tree$node.label)) {
      PCMTreeSetLabels(tree)
    } else if(length(unique(tree$node.label)) != length(tree$node.label)) {
      stop(
        paste(
          "PCMTree:: found duplicated labels in tree$node.label. This can cause",
          "likelihood calculation errors. ",
          "Please, ensure that all labels in node.label are unique or set ",
          "tree$node.label to NULL."))
    }

    N <- PCMTreeNumTips(tree)
    nodeLabels <- PCMTreeGetLabels(tree)

    if( !is.null(tree$edge.regime) ) {
      preorder <- PCMTreePreorder(tree)

      edge.part <- rep(1L, nrow(tree$edge))
      nextPart <- 2L
      nodesInOrder <- N + 1L

      # add an extra element in edge.regime denoting the regime for the root
      # as a one of its daughters' regimes picked at random
      tree$edge.regime <- c(
        tree$edge.regime,
        tree$edge.regime[preorder[1L]])

      endingAt <- order(rbind(tree$edge, c(0L, N + 1L))[, 2L])

      for(ei in preorder) {
        i <- tree$edge[ei, 2L]
        j <- tree$edge[ei, 1L]
        ej <- endingAt[j]

        if(tree$edge.regime[ei] != tree$edge.regime[ej]) {
          edge.part[ei] <- nextPart
          nextPart <- nextPart + 1L
          nodesInOrder <- c(nodesInOrder, i)
        } else if(j != N + 1L) {
          edge.part[ei] <- edge.part[ej]
        }
      }

      # now edge.part is integer vector containing the indices of the parts.
      # We replace these part-indices with their corresp. partition node-labels,
      # because these remain stable upon manipulation of the tree, such as insertion
      # of singleton nodes.
      edge.part <- nodeLabels[nodesInOrder[edge.part]]
      part.regime <- structure(
        tree$edge.regime[endingAt[nodesInOrder]],
        names = nodeLabels[nodesInOrder])

      tree$edge.part <- edge.part
      tree$part.regime <- part.regime

      clsTree <- class(tree)
      tree <- tree[-match("edge.regime", names(tree))]
      class(tree) <- clsTree

    } else {
      PCMTreeSetPartition(tree, N + 1L)
    }

    tree
  }
}

#' @importFrom ape as.phylo
#' @export
as.phylo.PCMTree <- function(x, ...) {
  x
}

#' Check that a tree is a PCMTree
#' @param tree a tree object.
#' @return a logical TRUE if `inherits(tree, "PCMTree")` is TRUE.
#'
#' @export
is.PCMTree <- function(tree) {
  inherits(tree, "PCMTree")
}

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

#' Set tip and internal node labels in a tree
#'
#' @param tree a phylo object or a PCMTree object. If this is a PCMTree object,
#' the internal edge.part and part.regime members will be updated accordingly.
#' @param labels a character vector in the order 1:PCMTreeNumNodes(tree) as
#' denoted in the tree$edge matrix.
#' @param inplace a logical indicating if the change should be done in place on
#' the object in the calling environment (in this case tree must not be a
#' temporary object, e.g. returned by another function call). Default is TRUE.
#'
#' @return if inplace is FALSE, a copy of tree with set or modified
#' tree$tip.label and tree$node.label. If the original tree has a member
#' edge.part, the returned tree has tree$edge.part and tree$part.regime updated.
#' If inplace is TRUE (the default), nothing is returned and the above changes
#' are made directly on the input tree.
#'
#' @seealso \code{\link{PCMTree}}
#'
#' @examples
#' tree <- ape::rtree(5)
#' tree$tip.label
#' # the following three are NULLs
#' tree$node.label
#' tree$edge.part
#' tree$part.regime
#'
#'
#' tree <- PCMTree(tree)
#' PCMTreeSetPartition(tree, c(6, 8))
#' tree$tip.label
#' tree$node.label
#' tree$edge.part
#' tree$part.regime
#'
#' PCMTreeSetLabels(
#'   tree, labels = paste0(c(rep("t", 5), rep("n", 4)), PCMTreeGetLabels(tree)))
#' PCMTreeGetLabels(tree)
#' tree$tip.label
#' tree$node.label
#' tree$edge.part
#' tree$part.regime
#'
#' @export
PCMTreeSetLabels <- function(
  tree, labels = as.character(1:PCMTreeNumNodes(tree)), inplace = TRUE) {

  if(!inherits(tree, "phylo")) {
    stop("PCMTreeSetLabels:: argument tree should be a phylo.")
  }

  N <- PCMTreeNumTips(tree)
  M <- PCMTreeNumNodes(tree)

  if(length(unique(labels)) != length(labels) ||
     length(labels) != M) {
    stop("PCMTreeSetLabels:: labels should contain M distinct elements." )
  }

  labels2 <- labels
  unname(labels2)
  if( !is.null(tree$node.label) ) {
    names(labels2) <- c(tree$tip.label, tree$node.label)
  } else {
    names(labels2) <- labels2
  }

  treeEdgePart <- tree$edge.part
  treePartRegime <- tree$part.regime

  if(inplace) {
    eval(substitute({
      if(!is.null(treeEdgePart)) {
          tree$edge.part <- unname(labels2[treeEdgePart])
        names(tree$part.regime) <- labels2[names(treePartRegime)]
      }
      tree$tip.label <- unname(labels2[seq_len(N)])
      tree$node.label <- unname(labels2[(N + 1):M])
    }), parent.frame())
  } else {
    if(!is.null(treeEdgePart)) {
      tree$edge.part <- unname(labels2[treeEdgePart])
      names(tree$part.regime) <- labels2[names(treePartRegime)]
    }
    tree$tip.label <- unname(labels2[seq_len(N)])
    tree$node.label <- unname(labels2[(N + 1):M])
    tree
  }
}

#' Get a vector of the tip and node labels in a tree
#' @param tree a phylo object
#' @return a character vector
#' @export
PCMTreeGetLabels <- function(tree) {
  if(!inherits(tree, "phylo")) {
    stop(
      "ERR:02650:PCMBase:PCMTree.R:PCMTreeGetLabels:: argument tree should be a phylo.")
  }
  if(is.null(tree$node.label)) {
    stop("ERR:02651:PCMBase:PCMTree.R:PCMTreeGetLabels:: the tree has no node.label member assigned.")
  }
  c(tree$tip.label, tree$node.label)
}

#' Get the node numbers associated with tip- or node-labels in a tree
#' @param tree a phylo object
#' @param labels a character vector with valid tip or node labels from tree
#' @return an integer vector giving the tip- or node- integer indices
#'  corresponding to labels.
#' @export
PCMTreeMatchLabels <- function(tree, labels) {
  allLabels <- PCMTreeGetLabels(tree)
  m <- match(labels, allLabels)
  if(any(is.na(m))) {
    stop("ERR:02660:PCMBase:PCMTree.R:PCMTreeMatchLabels:: some of the labels not found in PCMTreeGetLabels(tree).")
  }
  m
}


############################# Partition functions #############################

#' Number of unique parts on a tree
#' @param tree a phylo object
#' @return the number of different parts encountered on the tree branches
#' @export
PCMTreeNumParts <- function(tree) {
  tree <- PCMTree(tree)
  length(tree$part.regime)
}

#' @export
PCMRegimes.phylo <- function(obj) {
  PCMRegimes(PCMTree(obj))
}

#' @export
PCMRegimes.PCMTree <- function(obj) {
  unname(unique(obj$part.regime))
}

#' @export
PCMNumRegimes.phylo <- function(obj) {
  PCMNumRegimes(PCMTree(obj))
}

#' @export
PCMNumRegimes.PCMTree <- function(obj) {
  length(unique(obj$part.regime))
}

#' Get the starting branch' nodes for each part on a tree
#'
#' @details We call a starting branch the first branch from the root to the tips
#' of a given part. A starting node is the node at which a starting branch
#' ends.
#' @param tree a phylo object with an edge.part member denoting parts. The
#' function assumes that each part covers a linked set of branches on the tree.
#' @return a named integer vector with elements equal to the starting nodes for
#' each part in tree and names equal to the labels of these nodes.
#' @seealso \code{\link{PCMTreeSetPartition}}
#' @examples
#' set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
#' PCMTreeGetPartition(PCMTree(ape::rtree(20)))
#' @export
PCMTreeGetPartition <- function(tree) {
  tree <- PCMTree(tree)
  structure(
    match(names(tree$part.regime), PCMTreeGetLabels(tree)),
    names = names(tree$part.regime))
}

#' Set a partition of a tree by specifying the partition nodes
#'
#' @param tree a PCMTree object.
#' @param nodes a character vector containing tip or node labels or an integer
#' vector denoting tip or internal nodes in tree - the parts change at the
#' start of the branches leading to these nodes. Default:
#' c(PCMTreeNumTips(tree) + 1L).
#' @param inplace a logical indicating if the change should be done to the tree
#' in the calling environment (TRUE) or a copy of the tree with set edge.part
#' member should be returned (FALSE). Default is TRUE.
#' @return If inplace is TRUE nothing, otherwise a copy of the tree with set
#' edge.part member.
#'
#' @seealso \code{\link{PCMTreeGetPartition}}
#' @seealso \code{\link{PCMTree}}
#'
#' @examples
#' set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
#' tree <- PCMTree(ape::rtree(8))
#' PCMTreeSetLabels(tree, paste0("x", PCMTreeGetLabels(tree)))
#' PCMTreeGetPartition(tree)
#' PCMTreeGetPartNames(tree)
#' PCMTreeGetPartRegimes(tree)
#' \donttest{
#' PCMTreePlot(tree) + ggtree::geom_nodelab() + ggtree::geom_tiplab()
#' }
#' tree <- PCMTreeSetPartition(tree, c(12, 14), inplace = FALSE)
#' PCMTreeGetPartition(tree)
#' PCMTreeGetPartNames(tree)
#' PCMTreeGetPartRegimes(tree)
#' \donttest{
#' PCMTreePlot(tree) + ggtree::geom_nodelab() + ggtree::geom_tiplab()
#' }
#'
#' # reset the partition to a default one, where there is only one part:
#' PCMTreeSetPartition(tree)
#'
#' PCMTreeGetPartition(tree)
#' PCMTreeGetPartNames(tree)
#' PCMTreeGetPartRegimes(tree)
#' \donttest{
#' PCMTreePlot(tree) + ggtree::geom_nodelab() + ggtree::geom_tiplab()
#' }
#'
#' # reset the labels to the default labels which are character representations
#' # of the node ids
#' PCMTreeSetLabels(tree)
#' PCMTreeGetPartition(tree)
#' PCMTreeGetPartNames(tree)
#' PCMTreeGetPartRegimes(tree)
#'
#' @export
PCMTreeSetPartition <- function(
  tree, nodes = c(PCMTreeNumTips(tree) + 1L), inplace=TRUE) {

  if(!is.PCMTree(tree)) {
    stop(paste0(
      "PCMTreeSetPartition:: argument tree should be of class PCMTree."))
  }

  nodeLabels <- PCMTreeGetLabels(tree)

  if(is.character(nodes)) {
    nodes <- match(nodes, nodeLabels)
    if(any(is.na(nodes))) {
      stop(paste0(
        "PCMTreeSetPartition:: if nodes is a ",
        "character vector it should be a subset of PCMTreeGetLabels(tree),",
        "but some of the elements in nodes could not be matched."))
    }
  }

  preorder <- PCMTreePreorder(tree)

  edgePart <- rep(1L, nrow(tree$edge))
  nextPart <- 2L
  N <- PCMTreeNumTips(tree)
  nodesInOrder <- N + 1L

  endingAt <- order(rbind(tree$edge, c(0L, N + 1L))[, 2L])

  for(ei in preorder) {
    i <- tree$edge[ei, 2L]
    j <- tree$edge[ei, 1L]
    ej <- endingAt[j]

    if(i %in% nodes) {
      edgePart[ei] <- nextPart
      nextPart <- nextPart + 1L
      nodesInOrder <- c(nodesInOrder, i)
    } else if(j != N + 1L) {
      edgePart[ei] <- edgePart[ej]
    }
  }

  # now edgePart is integer vector containing the indices of the parts.
  # We replace these part-indices with their corresp. partition node-labels,
  # because these remain stable upon manipulation of the tree, such as
  # insertion of singleton nodes.
  edgePart <- nodeLabels[nodesInOrder[edgePart]]
  unname(edgePart)
  partRegime <- structure(
    seq_along(nodesInOrder), names = nodeLabels[nodesInOrder])

  if(inplace) {
    eval(substitute(tree$edge.part <- edgePart), parent.frame())
    eval(substitute(tree$part.regime <- partRegime), parent.frame())
  } else {
    tree$edge.part <- edgePart
    tree$part.regime <- partRegime
    tree
  }
}

#' Unique parts on a tree in the order of occurrence from the root to the tips (preorder)
#'
#' @param tree a phylo object with an additional member edge.part which should
#' be a character or an integer vector of length equal to the number of branches.
#' @return a character vector.
#' @export
PCMTreeGetPartNames <- function(tree) {
  tree <- PCMTree(tree)
  names(tree$part.regime)
}

#' Regime mapping for the parts in a tree
#' @param tree a PCMTree or a phylo object.
#' @return a named vector with names corresponding to the part names in tree
#' and values corresponding to regime names or ids.
#' @export
PCMTreeGetPartRegimes <- function(tree) {
  tree <- PCMTree(tree)
  tree$part.regime
}

#' Set regimes for the parts in a tree
#' @param tree a PCMTree object.
#' @param part.regime a named vector containing regimes to be assigned to
#' some of or to each of the parts in the tree.
#' @param setPartition a logical indicating if the partition of the tree should
#' be set as well. If this argument is set to TRUE, the names of part.regime
#' are passed as the nodes argument in a call to \code{PCMTreeSetPartition}.
#' Default: FALSE.
#' @param inplace a logical indicating if the change should be done to the tree
#' in the calling environment (TRUE) or a copy of the tree with set edge.part
#' member should be returned (FALSE). Default is TRUE.
#' @return If inplace is TRUE nothing, otherwise a copy of the tree with set
#' edge.part and part.regime members.
#'
#' @seealso \code{\link{PCMTree}}
#'
#' @examples
#' tree <- PCMTree(ape::rtree(25))
#' PCMTreeGetPartition(tree)
#' PCMTreeGetPartRegimes(tree)
#' PCMTreeGetPartNames(tree)
#'
#' PCMTreeSetPartRegimes(tree, c(`26` = 2))
#' PCMTreeGetPartition(tree)
#' PCMTreeGetPartRegimes(tree)
#' PCMTreeGetPartNames(tree)
#'
#' PCMTreeSetPartRegimes(tree, c(`26` = "global-regime"))
#' PCMTreeGetPartition(tree)
#' PCMTreeGetPartRegimes(tree)
#' PCMTreeGetPartNames(tree)
#'
#' # This should fail because no partition with nodes 26, 28 and 41 has been
#' # done.
#' ggplot2::should_stop(
#'   PCMTreeSetPartRegimes(tree, c(`26` = "a", `28` = "b", `41` = "c")))
#' # This should succeed and change the partition as well as regime assignment
#' PCMTreeSetPartRegimes(
#'   tree, c(`26` = "a", `28` = "b", `41` = "c"), setPartition = TRUE)
#' PCMTreeGetPartition(tree)
#' PCMTreeGetPartRegimes(tree)
#' PCMTreeGetPartNames(tree)
#'
#'
#'
#' set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
#' # number of tips
#' N <- 40
#'
#' # tree with one regime
#' tree.a <- ape::rtree(N)
#'
#' tree.a <- PCMTree(tree.a)
#'
#' PCMTreeSetPartRegimes(
#'   tree.a,
#'   part.regime = structure("a", names = as.character(N+1L)),
#'   setPartition = TRUE,
#'   inplace = TRUE)
#'
#'
#' \donttest{
#' PCMTreePlot(tree.a) + ggtree::geom_nodelab() + ggtree::geom_tiplab()
#' }
#' tree.ab <- tree.a
#' PCMTreeSetPartRegimes(
#'   tree.ab,
#'   part.regime = structure(c("a", "b"), names = as.character(c(N+1L, N+31L))),
#'   setPartition = TRUE,
#'   inplace = TRUE)
#' \donttest{
#' PCMTreePlot(tree.ab) + ggtree::geom_nodelab() + ggtree::geom_tiplab()
#' }
#' @export
PCMTreeSetPartRegimes <- function(
  tree, part.regime, setPartition = FALSE, inplace = TRUE) {

  if(!is.PCMTree(tree)) {
    stop(paste0(
      "PCMTreeSetPartRegimes:: argument tree should be of class PCMTree."))
  }

  if(is.null(names(part.regime)) ||
     length(names(part.regime)) != length(part.regime)) {
    stop("PCMTreeSetPartRegimes:: argument part.regime should be a named vector.")
  }

  tree2 <- tree

  if(setPartition) {
    tree2 <- PCMTreeSetPartition(tree2, names(part.regime), inplace = FALSE)
    treeEdgePart <- tree2$edge.part
  }
  if(any(is.na(match(names(part.regime), names(tree2$part.regime))))) {
    stop(paste0(
      "PCMTreeSetPartRegimes:: some of the names in the argument part.regime",
      "do not match with part names in the tree. "))
  }
  if(length(unique(names(part.regime))) < length(part.regime) ) {
    stop(paste0(
      "PCMTreeSetPartRegimes:: some part names in the argument part.regime are",
      "duplicated."))
  }

  treePartRegime <- tree2$part.regime
  treePartRegime[names(part.regime)] <- part.regime

  if(inplace) {
    if(setPartition) {
      eval( substitute(tree[["edge.part"]] <- treeEdgePart), parent.frame() )
    }
    eval( substitute(tree[["part.regime"]] <- treePartRegime), parent.frame() )
  } else {
    tree2$part.regime <- treePartRegime
    tree2
  }
}

#' Get the parts of the branches leading to a set of nodes or tips
#' @param tree a phylo object with an edge.part member denoting parts.
#' @param nodes an integer vector denoting the nodes.
#' Default is seq_len(PCMTreeNumNodes(tree).
#' @return a character vector denoting the parts of the branches
#' leading to the nodes, according to tree$edge.part.
#' @export
PCMTreeGetPartsForNodes <- function(
  tree, nodes = seq_len(PCMTreeNumNodes(tree))) {

  tree <- PCMTree(tree)
  parts <- tree$edge.part[match(nodes, tree$edge[, 2])]
  N <- PCMTreeNumTips(tree)
  parts[nodes == N+1L] <- tree$node.label[1]
  parts
}

#' Get the tips belonging to a part in a tree
#' @param tree a phylo object with an edge.regime member or a PCMTree object
#' @param part a character or integer denoting a part name in the tree.
#' @return an integer vector with the ids of the tips belonging to part
#' @examples
#' set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
#' tree <- ape::rtree(10)
#' regimes <- sample(letters[1:3], nrow(tree$edge), replace = TRUE)
#' PCMTreeSetRegimesForEdges(tree, regimes)
#' \donttest{
#' PCMTreePlot(tree) + ggtree::geom_nodelab() + ggtree::geom_tiplab()
#' }
#' part <- PCMTreeGetPartNames(tree)[1]
#' PCMTreeGetTipsInPart(tree, part)
#' print(part)
#'
#' @seealso \link{PCMTreeGetTipsInRegime}, \link{PCMTreeGetPartNames}, \link{PCMRegimes}, \link{PCMTreeGetPartRegimes}, \link{PCMTreeSetPartRegimes}
#' @export
PCMTreeGetTipsInPart <- function(tree, part) {

  tree <- PCMTree(tree)

  N <- PCMTreeNumTips(tree)
  tipEdges <- tree$edge[, 2] <= N & tree$edge.part == part
  tree$edge[tipEdges, 2]
}


#' Get the tips belonging to a regime in a tree
#' @param tree a phylo object with an edge.regime member or a PCMTree object
#' @param regime a character or integer denoting a regime in the tree.
#' @return an integer vector with the ids of the tips belonging to part
#' @examples
#' set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
#' tree <- ape::rtree(10)
#' regimes <- sample(letters[1:3], nrow(tree$edge), replace = TRUE)
#' PCMTreeSetRegimesForEdges(tree, regimes)
#' \donttest{
#' PCMTreePlot(tree) + ggtree::geom_nodelab() + ggtree::geom_tiplab()
#' }
#' regime <- PCMRegimes(tree)[1]
#' PCMTreeGetTipsInRegime(tree, regime)
#' print(regime)
#'
#' @seealso \link{PCMTreeGetTipsInPart}, \link{PCMTreeGetPartNames}, \link{PCMRegimes}, \link{PCMTreeGetPartRegimes}, \link{PCMTreeSetPartRegimes}, \link{PCMTreeGetPartition}
#'
#' @export
PCMTreeGetTipsInRegime <- function(tree, regime) {
  tipRegimes <- PCMTreeGetRegimesForNodes(tree)[seq_len(PCMTreeNumTips(tree))]
  which(regime == tipRegimes)
}


#' Model regimes (i.e. colors) associated with the branches in a tree
#' @param tree a PCMTree or a phylo object.
#' @return a vector with entries corresponding to the rows in tree$edge
#'   denoting the regime associated with each branch in the tree. The type
#'   of the vector element is defined by the type of tree$part.regime.
#' @export
PCMTreeGetRegimesForEdges <- function(tree) {
  tree <- PCMTree(tree)
  tree$part.regime[as.character(tree$edge.part)]
}

#' Set the regime for each individual edge in a tree explicitly
#' @note Calling this function overwrites the current partitioning of the tree.
#' @param tree a PCMTree or a phylo object.
#' @param regimes a vector of the length equal to `nrow(tree$edge)`.
#' @param inplace a logical indicating if the change should be done within the
#' tree in the calling environment or a copy of the tree with modified regime
#' assignment should be retrned.
#' @return if inplace is TRUE, nothing, otherwise a modified copy of tree.
#' @examples
#' set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
#' tree <- ape::rtree(10)
#' regimes <- sample(letters[1:3], nrow(tree$edge), replace = TRUE)
#' PCMTreeSetRegimesForEdges(tree, regimes)
#' \donttest{
#' PCMTreePlot(tree)
#' }
#'
#' @export
PCMTreeSetRegimesForEdges <- function(tree, regimes, inplace = TRUE) {
  if( !inherits(tree, "phylo") ) {
    stop(
      "PCMTreeSetRegimesForEdges: tree should be a PCMTree or a phylo object.")
  }

  if(!is.vector(regimes) || length(regimes) != nrow(tree$edge)) {
    stop(
      paste0(
        "PCMTreeSetRegimesForEdges: regimes should be a vector of length that ",
        "equal the number of rows in tree$edge."))
  }

  if(is.PCMTree(tree)) {
    # tree has edge.part and part.regime members which have to be overwritten.
    tree2 <- tree
    tree2$edge.part <- NULL
    tree2$part.regime <- NULL
    tree2 <- tree2[!sapply(tree2, is.null)]
    class(tree2) <- "phylo"
    tree2$edge.regime <- regimes

    tree2 <- PCMTree(tree2)

    if(inplace) {
      eval(substitute({
        tree$edge.part <- tree2$edge.part
          tree$part.regime <- tree2$part.regime
      }), parent.frame())
    } else {
      tree2
    }
  } else {
    if(inplace) {
      eval(substitute({
        tree$edge.regime <- regimes
      }), parent.frame())
    } else {
      tree$edge.regime <- regimes
      tree
    }
  }
}

#' Get the regimes of the branches leading to a set of nodes or tips
#' @param tree a phylo object with an edge.part member denoting parts.
#' @param nodes an integer vector denoting the nodes.
#' Default is seq_len(PCMTreeNumNodes(tree).
#' @return a character vector denoting the parts of the branches
#' leading to the nodes, according to tree$edge.part.
#' @importFrom data.table setkey
#' @export
PCMTreeGetRegimesForNodes <- function(
  tree, nodes = seq_len(PCMTreeNumNodes(tree) ) ) {

  # avoid nodes during check
  regime <- endNode <- NULL

  dtNodes <- PCMTreeDtNodes(tree)
  setkey(dtNodes, endNode)
  as.vector(dtNodes[list(nodes), regime])
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
#' @param skipNodesLabels a character vector of node labels in tree that should be
#' omitted as candidates for partition nodes. Default: character(0).
#' @param tableAncestors NULL (default) or an integer matrix returned by a
#' previous call to \code{PCMTreeTableAncestors(tree)}.
#' @param verbose a logical indicating if informative messages should be printed to
#' the console.
#'
#' @return a list of integer vectors.
#' @examples
#'
#' set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
#' tree <- PCMTree(ape::rtree(10))
#'
#' \donttest{
#' PCMTreePlot(tree) + ggtree::geom_nodelab() + ggtree::geom_tiplab()
#' }
#'
#' # list of all partitions into parts of at least 4 tips
#' PCMTreeListAllPartitions(tree, 4)
#'
#' # list of all partitions into parts of at least 3 tips
#' PCMTreeListAllPartitions(tree, 3)
#'
#' # list all partitions into parts of at least 3 tips, excluding the partitions
#' # where node 16 is one of the partition nodes:
#' PCMTreeListAllPartitions(tree, 3, "16")
#'
#' @export
PCMTreeListAllPartitions <- function(
  tree,
  minCladeSize,
  skipNodesLabels = character(),
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
    withoutNodesLabels = skipNodesLabels,
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

#' A data.table with time, part and regime information for the nodes in a tree
#' @param tree a phylo object with node-labels and parts
#' @return a data.table with a row for each node in tree and columns as follows:
#' \itemize{
#' \item{startNode }{the starting node of each edge or NA_integer_ for the root}
#' \item{endNode }{the end node of each edge or the root id for the root}
#' \item{startNodeLab }{the character label for the startNode}
#' \item{endNodeLab }{the character label for endNode}
#' \item{startTime }{the time (distance from the root node) for the startNode or
#'  NA_real_ for the root node}
#' \item{endTime }{the time (distance from the root node) for the endNode or
#'  NA_real_ for the root node}
#' \item{part }{the part to which the edge belongs, i.e. the part of the
#'  endNode}
#' \item{regime }{the regime to which the edge belongs, i.e. the regime of the
#'  part of the endNode}}
#' @importFrom data.table data.table rbindlist
#' @examples
#' PCMTreeDtNodes(PCMBaseTestObjects$tree.ab)
#' @export
PCMTreeDtNodes <- function(tree) {

  tree <- PCMTree(tree)

  N <- PCMTreeNumTips(tree)
  nodeTimes <- PCMTreeNodeTimes(tree)
  nodeLabels <- PCMTreeGetLabels(tree)

  dtForBranches <- data.table(
    startNode = tree$edge[, 1],
    endNode = tree$edge[, 2],
    startNodeLab = nodeLabels[tree$edge[, 1]],
    endNodeLab = nodeLabels[tree$edge[, 2]],
    startTime = nodeTimes[tree$edge[, 1]],
    endTime = nodeTimes[tree$edge[, 2]],
    part = tree$edge.part,
    regime = tree$part.regime[tree$edge.part])
  dtForRoot <- data.table(
    startNode = NA_integer_,
    endNode = PCMTreeNumTips(tree) + 1L,
    startNodeLab = NA_character_,
    endNodeLab = nodeLabels[PCMTreeNumTips(tree) + 1L],
    startTime = NA_real_,
    endTime = 0.0,
    part = PCMTreeGetPartsForNodes(tree, PCMTreeNumTips(tree) + 1L),
    regime = tree$part.regime[
      PCMTreeGetPartsForNodes(tree, PCMTreeNumTips(tree) + 1L)])
  rbindlist(list(dtForRoot, dtForBranches))
}

#' Prune the tree leaving one tip for each or some of its parts
#' @param tree a PCMTree or a phylo object.
#' @param partsToKeep a character vector denoting part names in the tree to be
#'  kept. Defaults to `PCMTreeGetPartNames(tree)`.
#'
#' @return a PCMTree object representing a pruned version of tree.
#'
#' @seealso PCMTreeSetPartition
#'
#' @importFrom data.table data.table
#' @importFrom ape drop.tip bind.tree
#'
#' @seealso PCMTree
#' @examples
#' set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
#' tree <- PCMTree(ape::rtree(25))
#' \donttest{
#' PCMTreePlot(tree) + ggtree::geom_nodelab(angle = 45) +
#'   ggtree::geom_tiplab(angle = 45)
#' }
#' backb <-  PCMTreeBackbonePartition(tree)
#' \donttest{
#' PCMTreePlot(backb) + ggtree::geom_nodelab(angle = 45) +
#'   ggtree::geom_tiplab(angle = 45)
#' }
#'
#' tree2 <- PCMTreeSetPartRegimes(
#'   tree, c(`26` = "a", `28` = "b"), setPartition = TRUE,
#'   inplace = FALSE)
#' \donttest{
#' PCMTreePlot(tree2) + ggtree::geom_nodelab(angle = 45) +
#'   ggtree::geom_tiplab(angle = 45)
#' }
#'
#' backb <-  PCMTreeBackbonePartition(tree2)
#' \donttest{
#' PCMTreePlot(backb) + ggtree::geom_nodelab(angle = 45) +
#'   ggtree::geom_tiplab(angle = 45)
#' }
#'
#' tree3 <- PCMTreeSetPartRegimes(
#'   tree, c(`26` = "a", `28` = "b", `41` = "c"), setPartition = TRUE,
#'   inplace = FALSE)
#' \donttest{
#' PCMTreePlot(tree3) + ggtree::geom_nodelab(angle = 45) +
#'   ggtree::geom_tiplab(angle = 45)
#' }
#' backb <-  PCMTreeBackbonePartition(tree3)
#' \donttest{
#' PCMTreePlot(backb) + ggtree::geom_nodelab(angle = 45) +
#'   ggtree::geom_tiplab(angle = 45)
#' }
#' backb41 <-  PCMTreeBackbonePartition(tree3, partsToKeep = "41")
#' \donttest{
#' PCMTreePlot(backb41) + ggtree::geom_nodelab(angle = 45) +
#'   ggtree::geom_tiplab(angle = 45)
#' }
#'
#' backbMoreNodes <- PCMTreeInsertSingletonsAtEpoch(
#'    backb, epoch = 3.7, minLength = 0.001)
#' PCMTreeGetPartRegimes(backbMoreNodes)
#' \donttest{
#' PCMTreePlot(backbMoreNodes) + ggtree::geom_nodelab(angle=45) +
#'   ggtree::geom_tiplab(angle=45)
#' }
#' backbMoreNodes <- PCMTreeInsertSingletonsAtEpoch(
#'    backbMoreNodes, epoch = 0.2, minLength = 0.001)
#' PCMTreeGetPartRegimes(backbMoreNodes)
#' \donttest{
#' PCMTreePlot(backbMoreNodes) + ggtree::geom_nodelab(angle=45) +
#'   ggtree::geom_tiplab(angle=45)
#' }
#' backbMoreNodes <- PCMTreeInsertSingletonsAtEpoch(
#'    backbMoreNodes, epoch = 1.2, minLength = 0.001)
#' PCMTreeGetPartRegimes(backbMoreNodes)
#'
#' \donttest{
#' PCMTreePlot(backbMoreNodes) + ggtree::geom_nodelab(angle=45) +
#'   ggtree::geom_tiplab(angle=45)
#' }
#' @export
PCMTreeBackbonePartition <- function(tree, partsToKeep = PCMTreeGetPartNames(tree)) {

  tree <- PCMTree(tree)

  # Needed to pass the R CMD CHECK.
  part <- endNode <- endTime <- NULL

  N <- PCMTreeNumTips(tree)
  nodeTimes <- PCMTreeNodeTimes(tree)
  nodeLabels <- PCMTreeGetLabels(tree)

  partNodesTree <- PCMTreeGetPartition(tree)

  dtBranchParts <- data.table(
    startNode = tree$edge[, 1], endNode = tree$edge[, 2],
    startTime = nodeTimes[tree$edge[, 1]], endTime = nodeTimes[tree$edge[, 2]],
    part = tree$edge.part)

  dtTipParts <- dtBranchParts[endNode <= N]
  dtTipParts <- dtTipParts[
    , list(endNode=endNode[which.max(endTime)]),
    keyby = part][list(partsToKeep)]

  tipsToDrop <- setdiff(seq_len(N), dtTipParts[, endNode])

  # attach a dummy tip to the root of the tree to prevent dropping the root in
  # some cases (e.g. when there is only one part in the tree).
  dummyTipTree <- list(edge = matrix(c(2, 1), nrow = 1L),
                       tip.label = "----dymmy---tip----",
                       edge.length = 1.0,
                       Nnode = 1L)

  class(dummyTipTree) <- "phylo"
  tree <- bind.tree(tree, dummyTipTree, N + 1L)

  tree2 <- drop.tip(tree, tipsToDrop, collapse.singles = FALSE)

  # now, remove the dummy tip from tree2:
  dummyTipId <- match(dummyTipTree$tip.label, tree2$tip.label)
  tree2$edge.length <- tree2$edge.length[-match(dummyTipId, tree2$edge[,2])]
  tree2$edge <- tree2$edge[-match(dummyTipId, tree2$edge[,2]), , drop = FALSE]
  tree2$edge <- apply(tree2$edge, 1:2, function(n) if(n > dummyTipId) n-1 else n)
  tree2$tip.label <- tree2$tip.label[-dummyTipId]

  tree2 <- PCMTree(tree2)

  nodeLabels2 <- PCMTreeGetLabels(tree2)
  partNodesTree2 <- match(nodeLabels[partNodesTree], nodeLabels2)
  partNodesTree2 <- partNodesTree2[!is.na(partNodesTree2)]
  PCMTreeSetPartRegimes(
    tree = tree2,
    part.regime = structure(
      tree$part.regime[nodeLabels2[partNodesTree2]],
      names = nodeLabels2[partNodesTree2]),
    setPartition = TRUE,
    inplace = TRUE)
  tree2
}


#' Which couples from a given set of nodes in a tree belong to the same part?
#'
#' @param tree a PCMTree object or a phylo object.
#' @param nodes an integer vector of length L >= 2 denoting a set of nodes in
#'  the tree.
#' @param upperTriangle logical indicating if all duplicated entries and
#'  diagonal entries sould be set to NA (by default TRUE).
#' @param returnVector logical indicating if a vector instead of a matrix
#'  should be returned (corresponding to calling as.vector on the resulting
#'  matrix and removing
#' NAs). Default: TRUE
#' @return a L x L logical matrix with TRUE on the diagonal and for each couple
#' of tips that belong to the same part or regime. If returnVector is TRUE
#' (default) only a vector of the non-NA entries will be returned.
#' @examples
#' set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
#' tree <- PCMTree(ape::rtree(8))
#' PCMTreeMatrixNodesInSamePart(tree, returnVector = FALSE)
#'
#' PCMTreeSetPartition(tree, c(10, 12))
#' PCMTreeMatrixNodesInSamePart(tree, returnVector = FALSE)
#'
#' PCMTreeMatrixNodesInSamePart(tree)
#' PCMTreeMatrixNodesInSamePart(tree, seq_len(PCMTreeNumTips(tree)))
#' PCMTreeMatrixNodesInSamePart(
#'   tree, seq_len(PCMTreeNumTips(tree)), returnVector = FALSE)
#'
#' @export
PCMTreeMatrixNodesInSamePart <- function(
  tree, nodes = seq_len(PCMTreeNumNodes(tree)),
  upperTriangle = TRUE, returnVector = TRUE) {

  tree <- PCMTree(tree)

  nodeParts <- PCMTreeGetPartsForNodes(tree, nodes)

  L <- length(nodeParts)

  mat <- matrix(FALSE, length(nodes), length(nodes))
  colnames(mat) <- rownames(mat) <- PCMTreeGetLabels(tree)[nodes]
  diag(mat) <- TRUE

  for(i in seq_len(L)) {
    for(j in seq_len(i - 1)) {
      mat[i,j] <- mat[j,i] <- nodeParts[i] == nodeParts[j]
    }
  }

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

#' Which couples from a given set of nodes in a tree belong to the same regime?
#'
#' @rdname PCMTreeMatrixNodesInSamePart
#'
#' @return a L x L logical matrix with TRUE on the diagonal and for each couple
#' of tips that belong to the same part or regime. If returnVector is TRUE (default)
#' only a vector of the non-NA entries will be returned.
#' @examples
#' set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
#' tree <- PCMTree(ape::rtree(8))
#' PCMTreeMatrixNodesInSamePart(tree, returnVector = FALSE)
#'
#' PCMTreeSetPartition(tree, c(10, 12))
#' PCMTreeMatrixNodesInSamePart(tree, returnVector = FALSE)
#'
#' PCMTreeMatrixNodesInSamePart(tree)
#' PCMTreeMatrixNodesInSamePart(tree, seq_len(PCMTreeNumTips(tree)))
#' PCMTreeMatrixNodesInSamePart(
#'   tree, seq_len(PCMTreeNumTips(tree)), returnVector = FALSE)
#'
#' @export
PCMTreeMatrixNodesInSameRegime <- function(
  tree, nodes = seq_len(PCMTreeNumNodes(tree)),
  upperTriangle = TRUE, returnVector = TRUE) {

  tree <- PCMTree(tree)

  nodeRegimes <- PCMTreeGetPartRegimes(tree)[
    PCMTreeGetPartsForNodes(tree, nodes)]

  L <- length(nodeRegimes)

  mat <- matrix(FALSE, length(nodes), length(nodes))
  colnames(mat) <- rownames(mat) <- PCMTreeGetLabels(tree)[nodes]
  diag(mat) <- TRUE

  for(i in seq_len(L)) {
    for(j in seq_len(i - 1)) {
      mat[i,j] <- mat[j,i] <- nodeRegimes[i] == nodeRegimes[j]
    }
  }

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

#' Pre-order tree traversal
#' @param tree a phylo object with possible singleton nodes (i.e. internal
#'  nodes with one daughter node)
#' @return a vector of indices of edges in tree$edge in pre-order.
#' @export
PCMTreePreorder <- function(tree) {
  # number of tips
  N <- length(tree$tip.label)

  # total number of nodes in the tree is the number of edges + 1 for the root
  M <- dim(tree$edge)[1]+1

  ordFrom <- order(tree$edge[,1])

  # we need the ordered edges in order to easily traverse all edges starting
  # from a given node
  iFrom <- match(seq_len(M), tree$edge[ordFrom, 1])

  # the result is a vector of edge indices in the breadth-first search order
  res <- vector(mode='integer', length = M-1)

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
#' @param tree a phylo object with possible singleton nodes (i.e. internal nodes
#' with one daughter node)
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
#' @param preorder an integer vector returned by a previous call to
#' \code{PCMTreePreorder(tree)}. Default \code{PCMTreePreorder(tree)}.
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
#' @return A vector of size the number of nodes in the tree (tips, root,
#'   internal) containing the time from the root to the corresponding node in
#'   the tree.
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

#' Slit a tree at a given internal node into a clade rooted at this node and the remaining tree after dropping this clade
#' @param tree a PCMTree object.
#' @param node an integer or character indicating a root, internal or tip node
#' @param tableAncestors an integer matrix returned by a previous call to
#'  PCMTreeTableAncestors(tree) or NULL.
#' @param X an optional k x N matrix with trait value vectors for each tip in
#'  tree.
#' @return A list containing two named phylo objects:
#' \itemize{
#' \item{clade }{The subtree (clade) starting at \code{node}.}
#' \item{Xclade }{The portion of X attributable to the tips in clade; NULL if X is NULL.}
#' \item{rest }{The tree resulting after dropping all tips in the clade.}
#' \item{Xrest }{The portion of X attributable to the tips in rest; NULL if X is NULL.}
#' }
#' @details In the current implementation, the edge.jump and edge.part members
#' of the tree will be discarded and not present in the clade.
#'
#' @importFrom ape drop.tip
#' @examples
#' set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
#' tree <- PCMTree(ape::rtree(25))
#' \donttest{
#' PCMTreePlot(tree) + ggtree::geom_nodelab(angle = 45) +
#'   ggtree::geom_tiplab(angle = 45)
#' }
#' spl <- PCMTreeSplitAtNode(tree, 28)
#' \donttest{
#' PCMTreePlot(PCMTree(spl$clade)) + ggtree::geom_nodelab(angle = 45) +
#'   ggtree::geom_tiplab(angle = 45)
#' }
#' \donttest{
#' PCMTreePlot(PCMTree(spl$rest)) + ggtree::geom_nodelab(angle = 45) +
#'   ggtree::geom_tiplab(angle = 45)
#' }
#' @export
PCMTreeSplitAtNode <- function(tree, node, tableAncestors = PCMTreeTableAncestors(tree), X=NULL) {

  if(!inherits(tree, "phylo")) {
    stop("PCMTreeSplit:: tree must be a PCMTree object.")
  }

  N <- PCMTreeNumTips(tree)
  M <- PCMTreeNumNodes(tree)
  nodeLabels <- PCMTreeGetLabels(tree)


  if(is.character(node)) {
    nodeLab <- node
    node <- match(nodeLab, nodeLabels)
    if(is.na(node)) {
      stop(paste0(
        "ERR:02601:PCMBase:PCMTree.R:PCMTreeSplit:: character node (",
        node, ") was not matched against the nodeLabels in tree."))
    }
  } else {
    nodeLab <- try(nodeLabels[node], silent = TRUE)
    if(class(nodeLab) == "try-error") {
      stop(paste0(
        "ERR:02602:PCMBase:PCMTree.R:PCMTreeSplit:: non-character type node (",
        node, ") was not found: ", nodeLab, "."))
    }
  }

  if(node == N+1) {
    list(clade = tree,
         Xclade = X,
         rest = NULL,
         Xrest = NULL)
  } else {
    partRegimes <- PCMTreeGetPartRegimes(tree)
    regimeSplitNode <- partRegimes[PCMTreeGetPartsForNodes(tree, node)]
    names(regimeSplitNode) <- nodeLab

    # remove edge.part and part.regime (or in old tree objects, edge.regime),
    # and edge.jump from the tree - we currently do not update them.
    if( !is.null(tree$edge.regime) ) {
      tree <- tree[- which(names(tree) == "edge.regime")]
      class(tree) <- "phylo"
    }
    if( !is.null(tree$edge.part) ) {
      tree <- tree[- which(names(tree) == "edge.part")]
      class(tree) <- "phylo"
    }
    if( !is.null(tree$part.regime) ) {
      tree <- tree[- which(names(tree) == "part.regime")]
      class(tree) <- "phylo"
    }
    if(!is.null(tree$edge.jump)) {
      tree <- tree[- which(names(tree) == "edge.jump")]
      class(tree) <- "phylo"
    }

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
      colnames(X) <- tree$tip.label <- as.character(seq_len(N))
    }

    clade = drop.tip(
      tree,
      tip = setdiff(seq_len(N), tipsClade),
      trim.internal = TRUE, collapse.singles = FALSE)

    # Because we used collapse.singles=FALSE, clade is still holding the entire
    # lineage from the root to node. We need to cut that part manually so that
    # clade is rooted at node. We don't have to do this for rest, because it
    # should indeed be rooted at the tree-root.
    cladeNodeLabels <- PCMTreeGetLabels(clade)
    cladeEdgeNodeLabels <-
      cbind(cladeNodeLabels[clade$edge[, 1]], cladeNodeLabels[clade$edge[, 2]])
    cladeRootLab <- cladeNodeLabels[PCMTreeNumTips(clade) + 1]

    # TODO:  avoid this while loop and remove all root-nodes at once
    while(cladeRootLab != nodeLab) {
      ei <- which(cladeEdgeNodeLabels[, 1] == cladeRootLab)

      cladeRootLabNew <- cladeEdgeNodeLabels[ei, 2]
      clade$edge <- clade$edge[-ei, ]
      clade$edge[clade$edge > length(tipsClade)] <-
        clade$edge[clade$edge > length(tipsClade)] - 1
      clade$edge.length <- clade$edge.length[-ei]
      cladeEdgeNodeLabels <- cladeEdgeNodeLabels[-ei, ]
      clade$node.label <-
        clade$node.label[-match(cladeRootLab, clade$node.label)]
      clade$Nnode <- clade$Nnode - 1
      cladeRootLab <- cladeRootLabNew
    }

    rest = drop.tip(
      tree,
      tip = tipsClade,
      trim.internal = TRUE,
      collapse.singles = FALSE)

    # restore the partition and regimes in clade and rest
    clade <- PCMTree(clade)
    partRegimesClade <-
      partRegimes[intersect(names(partRegimes), PCMTreeGetLabels(clade))]
    if( !(nodeLab %in% names(partRegimesClade)) ) {
      partRegimesClade <- c(regimeSplitNode, partRegimesClade)
    }
    PCMTreeSetPartRegimes(clade, partRegimesClade, setPartition = TRUE)

    rest <- PCMTree(rest)
    partRegimesRest <-
      partRegimes[intersect(names(partRegimes), PCMTreeGetLabels(rest))]
    PCMTreeSetPartRegimes(rest, partRegimesRest, setPartition = TRUE)

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
#' @param tree a PCMTree object.
#' @param cladeRootNode a character string denoting the label or an integer denoting a node in the tree.
#' @param tableAncestors an integer matrix returned by a previous call to PCMTreeTableAncestors(tree) or NULL.
#' @param X an optional k x N matrix with trait value vectors for each tip in tree.
#' @param returnList logical indicating if only the phylo object associated
#'  with the clade should be returned. Defaults to \code{!is.null(X)}
#' @return If returnList is FALSE, a phylo object associated with the clade,
#'  otherise, a list with two named members :
#' \itemize{
#' \item{tree}{the phylo object associated with the clade}
#' \item{X}{the submatrix of X with columns corresponding to the tips in the clade}
#' }
#' @seealso PCMTreeSpliAtNode PCMTreeDropClade
#' @examples
#' set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
#' tree <- PCMTree(ape::rtree(25))
#' PCMTreeSetPartRegimes(
#'   tree, c(`26`="a", `28`="b", `45`="c"), setPartition = TRUE)
#' \donttest{
#' PCMTreePlot(tree, palette=c(a = "red", b = "green", c = "blue")) +
#'   ggtree::geom_nodelab(angle = 45) + ggtree::geom_tiplab(angle = 45)
#' }
#' blueTree <- PCMTreeExtractClade(tree, 45)
#' PCMTreeGetPartRegimes(blueTree)
#' \donttest{
#' PCMTreePlot(blueTree, palette=c(a = "red", b = "green", c = "blue")) +
#'   ggtree::geom_nodelab(angle = 45) + ggtree::geom_tiplab(angle = 45)
#' }
#' # we need to use the label here, because the node 29 in tree is not the same
#' # id in redGreenTree:
#' blueTree2 <- PCMTreeDropClade(blueTree, "48")
#' \donttest{
#' PCMTreePlot(blueTree2, palette=c(a = "red", b = "green", c = "blue")) +
#'   ggtree::geom_nodelab(angle = 45) + ggtree::geom_tiplab(angle = 45)
#' }
#'
#' @export
PCMTreeExtractClade <- function(tree, cladeRootNode, tableAncestors = NULL, X=NULL, returnList = !is.null(X)) {
  if(is.character(cladeRootNode)) {
    if(!is.character(tree$node.label)) {
      stop(paste0(
        "PCMTreeExtractClade::", cladeRootNode,
        ": cladeRootNode is a character string but tree$node.label is missing or not a character vector."))
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
        stop(paste0(
          "PCMTreeExtractClade::", cladeRootNode,
          ": cladeRootNode of character-type was not found in tree$node.label"))
      }
    }
  } else {
    cladeRootNodeNumber <- as.integer(cladeRootNode)
    if(cladeRootNodeNumber <= 0 || cladeRootNodeNumber > PCMTreeNumNodes(tree)) {
      stop(paste0(
        "PCMTreeExtractClade::", cladeRootNode,
        ": cladeRootNode of integer should be between 1 and M=",
        PCMTreeNumNodes(tree), " (the number of nodes in the tree)."))
    }
  }

  spl <- PCMTreeSplitAtNode(
    tree, cladeRootNodeNumber, tableAncestors = tableAncestors, X = X)

  if(!returnList) {
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
#' @param returnList logical indicating if a list of the phylo object
#' associated with the tree after dropping the clade and the corresponding
#' entries in X should be returned. Defaults to \code{!is.null(X)}
#' @param errorOnMissing logical indicating if an error should be rased if
#' cladeRootNode is not among the nodes in tree. Default FALSE, meaning that if
#' cladeRootNode is not a node in tree the tree (and X if
#' returnList is TRUE) is/are returned unchanged.
#' @return If returnList is FALSE, a phylo object associated with the remaining
#' tree after dropping the clade, otherise, a list with two named members :
#' \itemize{
#' \item{tree}{the phylo object associated with the remaining tree after dropping the clade}
#' \item{X}{the submatrix of X with columns corresponding to the tips in the remaining tree}
#' }
#' @seealso PCMTreeSpliAtNode PCMTreeExtractClade
#' @examples
#' set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
#' tree <- PCMTree(ape::rtree(25))
#' PCMTreeSetPartRegimes(
#'   tree, c(`26`="a", `28`="b", `45`="c"), setPartition = TRUE)
#' \donttest{
#' PCMTreePlot(tree, palette=c(a = "red", b = "green", c = "blue")) +
#'   ggtree::geom_nodelab(angle = 45) + ggtree::geom_tiplab(angle = 45)
#' }
#' redGreenTree <- PCMTreeDropClade(tree, 45)
#' PCMTreeGetPartRegimes(redGreenTree)
#' \donttest{
#' PCMTreePlot(redGreenTree, palette=c(a = "red", b = "green", c = "blue")) +
#'   ggtree::geom_nodelab(angle = 45) + ggtree::geom_tiplab(angle = 45)
#' }
#' # we need to use the label here, because the node 29 in tree is not the same
#' # id in redGreenTree:
#' redGreenTree2 <- PCMTreeDropClade(redGreenTree, "29")
#' \donttest{
#' PCMTreePlot(redGreenTree2, palette=c(a = "red", b = "green", c = "blue")) +
#'   ggtree::geom_nodelab(angle = 45) + ggtree::geom_tiplab(angle = 45)
#' }
#'
#' @export
PCMTreeDropClade <- function(tree, cladeRootNode, tableAncestors = NULL, X=NULL, returnList = !is.null(X), errorOnMissing = FALSE) {
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

  res <- if(!returnList) {
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
    spl <- PCMTreeSplitAtNode(
      tree, cladeRootNodeNumber, tableAncestors = tableAncestors, X = X)

    if(!returnList) {
      res <- spl$rest
    } else {
      res <- list(tree = spl$rest, X = spl$Xrest)
    }
  }

  res
}


#' Perfrorm nested extractions or drops of clades from a tree
#' @param tree a phylo object with named tips and internal nodes
#' @param expr a character string representing an R expression of nested calls
#'  of functions
#' \code{E(x,node)} denoting extracting the clade rooted at node from the tree
#'  x, or \code{D(x,node)}, denoting dropping the clade rooted at node from the
#'  tree x. These calls can be nested, i.e. x can be either the symbol x
#'  (corresponding to the original tree passed as argument) or a nested call to
#'  d or e.
#' @return the resulting phylo object from evaluating expr on tree.
#'
#' @examples
#' set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
#' tree <- PCMTree(ape::rtree(25))
#' PCMTreeSetPartRegimes(
#'   tree, c(`26`="a", `28`="b", `45`="c", `47`="d"), setPartition = TRUE)
#' \donttest{
#' PCMTreePlot(
#'   tree, palette=c(a = "red", b = "green", c = "blue", d = "magenta")) +
#'   ggtree::geom_nodelab(angle = 45) + ggtree::geom_tiplab(angle = 45)
#' }
#' bluePart <- PCMTreeEvalNestedEDxOnTree("D(E(tree,45),47)", tree)
#' PCMTreeGetPartRegimes(bluePart)
#' \donttest{
#' PCMTreePlot(
#'   bluePart, palette=c(a = "red", b = "green", c = "blue", d = "magenta")) +
#'   ggtree::geom_nodelab(angle = 45) + ggtree::geom_tiplab(angle = 45)
#' }
#'
#' # Swapping the D and E calls has the same result:
#' bluePart2 <- PCMTreeEvalNestedEDxOnTree("E(D(tree,47),45)", tree)
#' PCMTreeGetPartRegimes(bluePart2)
#' \donttest{
#' PCMTreePlot(
#'   bluePart2, palette=c(a = "red", b = "green", c = "blue", d = "magenta")) +
#'   ggtree::geom_nodelab(angle = 45) + ggtree::geom_tiplab(angle = 45)
#' }
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
#' @param threshold a positive numeric; only branches longer than threshold
#'  will be returned; Default 0.
#' @return a named list with an integer vector element "nodes" denoting the
#' ending nodes for each branch crossing epoch and numeric vector element
#' "positions" denoting the root-ward offset from each node in nodes.
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
#' @param nodes an integer vector denoting the terminating nodes of the edges
#'  on which a singleton node is to be inserted. This vector should not have
#'  duplicated nodes - if there is a need to insert two or more singleton nodes
#'  at distinct positions of the same edge, this should be done by calling the
#'  function several times with the longest position first and so on .
#' @param positions a positive numeric vector of the same length as nodes
#'  denoting the root-ward distances from nodes at which the singleton nodes
#'  should be inserted.
#' @param epoch a numeric indicating a distance from the root at which a
#' singleton node should be inserted in all lineages that are alive at that
#' time.
#' @param minLength a numeric indicating the minimum allowed branch-length
#' after dividing a branch by insertion of a singleton nodes. No singleton node
#' is inserted if this would result in a branch shorter than `minLength`. Note
#' that this condition is checked only in `PCMTreeInsertSingletonsAtEpoch`.
#'
#' @importFrom ape bind.tree drop.tip
#' @return a modified version of tree with inserted singleton nodes at the
#'  specified locations
#' @seealso \code{\link{PCMTreeEdgeTimes}} \code{\link{PCMTreeLocateEpochOnBranches}} \code{\link{PCMTreeLocateMidpointsOnBranches}}
#' @examples
#' set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
#' tree <- PCMTree(ape::rtree(25))
#' PCMTreeSetPartRegimes(
#'   tree, c(`26`="a", `28`="b", `45`="c", `47`="d"), setPartition = TRUE)
#' \donttest{
#' PCMTreePlot(
#'   tree,
#'   palette=c(a = "red", b = "green", c = "blue", d = "magenta")) +
#'   ggtree::geom_nodelab(angle = 45) + ggtree::geom_tiplab(angle = 45)
#' }
#' cbind(tree$edge, PCMTreeEdgeTimes(tree))
#'
#' id47 <- PCMTreeMatchLabels(tree, "47")
#' length47 <- PCMTreeGetBranchLength(tree, id47)
#'
#' # insert a singleton at 0.55 (root-ward) from node 47
#' tree <- PCMTreeInsertSingletons(tree, nodes = "47", positions = (length47/2))
#' \donttest{
#' PCMTreePlot(
#'   tree,
#'   palette=c(a = "red", b = "green", c = "blue", d = "magenta")) +
#'   ggtree::geom_nodelab(angle = 45) + ggtree::geom_tiplab(angle = 45)
#' }
#' # this fails, because the branch leading to node "47" is shorter now (0.55).
#' ggplot2::should_stop(
#'   tree <- PCMTreeInsertSingletons(
#'     tree, nodes = "47", positions= 2* length47 / 3))
#'
#' # the tree is the same
#' \donttest{
#' PCMTreePlot(
#'   tree, palette=c(a = "red", b = "green", c = "blue", d = "magenta")) +
#'   ggtree::geom_nodelab(angle = 45) + ggtree::geom_tiplab(angle = 45)
#' }
#' # we can insert at a position within the edge:
#' tree <- PCMTreeInsertSingletons(tree, nodes = "47", positions = length47/3)
#' \donttest{
#' PCMTreePlot(
#'   tree, palette=c(a = "red", b = "green", c = "blue", d = "magenta")) +
#'   ggtree::geom_nodelab(angle = 45) + ggtree::geom_tiplab(angle = 45)
#' }
#'
#' # Insert singletons at all branches crossing a given epoch. This will skip
#' # inserting singleton nodes where the resulting branches would be shorter
#' # than 0.1.
#' tree <- PCMTreeInsertSingletonsAtEpoch(tree, 2.3)
#' \donttest{
#' PCMTreePlot(
#'   tree, palette=c(a = "red", b = "green", c = "blue", d = "magenta")) +
#'   ggtree::geom_nodelab(angle = 45) + ggtree::geom_tiplab(angle = 45)
#' }
#' # Insert singletons at all branches crossing a given epoch
#' tree <- PCMTreeInsertSingletonsAtEpoch(tree, 2.3, minLength = 0.001)
#' \donttest{
#' PCMTreePlot(
#'   tree,
#'   palette=c(a = "red", b = "green", c = "blue", d = "magenta")) +
#'   ggtree::geom_nodelab(angle = 45) + ggtree::geom_tiplab(angle = 45)
#' }
#' @export
PCMTreeInsertSingletons <- function(tree, nodes, positions) {

  tree <- PCMTree(tree)

  M <- PCMTreeNumNodes(tree)
  N <- PCMTreeNumTips(tree)
  originalLabels <- PCMTreeGetLabels(tree)

  if(is.character(nodes)) {
    nodes <- match(nodes, originalLabels)
  }
  if(any(is.na(nodes)) ||
     max(nodes, na.rm = TRUE) > PCMTreeNumNodes(tree) ||
     min(nodes, na.rm = TRUE) < 1L) {
    stop(paste0(
      "PCMTreeInsertSingletons:: Some of the nodes were NA or they could not",
      "be matched against nodes in tree."))
  }

  PCMTreeSetLabels(tree, paste0("x_", PCMTreeGetLabels(tree)))

  names(originalLabels) <- PCMTreeGetLabels(tree)
  partRegimeWithXLabels <- PCMTreeGetPartRegimes(tree)

  edgeNames <- PCMTreeGetLabels(tree)[tree$edge[, 2]]

  edgeTimes <- PCMTreeEdgeTimes(tree)
  rownames(edgeTimes) <- edgeNames

  # which edges should be processed
  names(nodes) <- names(positions) <- PCMTreeGetLabels(tree)[nodes]
  edgesToCut <- tree$edge[, 2] %in% nodes
  edgesToCutNames <- edgeNames[edgesToCut]

  tipTree <- list(edge = matrix(c(2, 1), nrow = 1, ncol = 2),
                  tip.label = "y_1",
                  node.label = "y_2",
                  edge.length = c(1.0),
                  Nnode = 1)
  class(tipTree) <- "phylo"

  for(edgeName in edgesToCutNames) {
    # find the number of the ending node for the edge. DON'T USE CASHED LABELS
    # HERE, but use PCMTreeGetLabels(tree)!
    node <- match(edgeName, PCMTreeGetLabels(tree))

    # find the position (root-ward offset from node) where to cut

    # add, then drop the tipTree using ape's functions
    tree <- bind.tree(tree, tipTree, node, positions[edgeName])
    tree <- drop.tip(tree, "y_1", collapse.singles = FALSE)

    # inserted edge
    edgeNameNew <- paste0(
      "i_",
      round(edgeTimes[edgeName, 2] - positions[edgeName], 2), "_", edgeName)

    tree$node.label[is.na(tree$node.label)] <- edgeNameNew

    # The edge leading to the new internal node take the same part as ITS
    # DAUGHTER edge.
    # If the edge we just cut was leading to a partition node, then this
    # node is no more a partition node, because the newly inserted node is
    # its parent and it has to becomes the partition node.
    # Also the part should be renamed. The regime for the (now renamed) part is
    # preserved.
    idxPartNodeInPartRegimes <- match(edgeName, names(partRegimeWithXLabels))
    if(!is.na(idxPartNodeInPartRegimes)) {
      names(partRegimeWithXLabels)[idxPartNodeInPartRegimes] <- edgeNameNew
    }
  }

  PCMTreeSetPartRegimes(
    tree, part.regime = partRegimeWithXLabels, setPartition = TRUE)

  # restore original node labels
  restoredOriginalLabels <- PCMTreeGetLabels(tree)
  m <- match(restoredOriginalLabels, names(originalLabels))
  restoredOriginalLabels[!is.na(m)] <- originalLabels[m[!is.na(m)]]

  # restore original node-names for all old nodes
  PCMTreeSetLabels(tree, labels = unname(restoredOriginalLabels))

  tree
}

#' @describeIn PCMTreeInsertSingletons
#'
#' @export
PCMTreeInsertSingletonsAtEpoch <- function(tree, epoch, minLength = 0.1) {

  tree <- PCMTree(tree)

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
#' @param epoch a non-negative numeric
#' @return an integer vector
#' @export
PCMTreeNearestNodesToEpoch <- function(tree, epoch) {

  if(epoch == 0) {
    ans <- PCMTreeNumTips(tree) + 1L
  } else {
    nodeTimes <- PCMTreeNodeTimes(tree)
    if(epoch > max(nodeTimes)) {
      epoch <- max(nodeTimes)
    }
    points <- PCMTreeLocateEpochOnBranches(tree, epoch)
    ans <- sapply(seq_along(points$nodes), function(i) {
      node <- points$nodes[i]
      parentNode <- tree$edge[tree$edge[, 2] == points$nodes[i], 1]

      if(points$positions[i] > (nodeTimes[points$nodes[i]] - points$positions[i] - nodeTimes[parentNode])) {
        parentNode
      } else {
        node
      }
    })
  }
  unique(ans)
}

#' A character representation of a phylo object.
#'
#' @param tree a phylo object.
#' @param includeLengths logical. Default: FALSE.
#' @param includePartition logical. Default: FALSE.
#' @return a character string.
#' @export
PCMTreeToString <- function(
  tree, includeLengths = FALSE, includePartition = FALSE) {

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
  if(includePartition) {
    startingNodesParts <- PCMTreeGetPartition(tree)
    startingNodesPartsLabels <- nodeLabs[startingNodesParts]
    attributes(startingNodesPartsLabels) <- NULL
    unname(startingNodesPartsLabels)
  } else {
    startingNodesPartsLabels <- ""
  }

  paste0(toString(edgeLabelsOrdered), "; ",
         toString(edgeLengthsOrdered), "; ",
         toString(startingNodesPartsLabels))
}



#' Plot a tree with parts and regimes assigned to these parts
#'
#' @param tree a PCMTree or a phylo object.
#' @param palette a named vector of colors corresponding to the regimes in tree
#' @param ... Arguments passed to ggtree, e.g. layout = 'fan', open.angle = 8,
#' size=.25.
#' @note This function requires that the ggtree package is installed.
#' At the time of releasing this version the ggtree package is not available on
#' CRAN. Check the
#' \href{https://guangchuangyu.github.io/software/ggtree/}{ggtree homepage} for
#' instruction on how to install this package:
#' .
#' @importFrom data.table data.table
#' @importFrom grDevices hcl
#' @importFrom ggplot2 aes scale_color_manual
#'
#' @export
PCMTreePlot <- function(
  tree,
  palette = PCMColorPalette(PCMNumRegimes(tree),
                            PCMRegimes(tree)), ...) {

  if(!is.PCMTree(tree)) {
    tree <- PCMTree(tree)
  }

  # Needed to pass the check
  regime <- NULL

  if(requireNamespace("ggtree")) {
    N <- PCMTreeNumTips(tree)
    R <- PCMTreeNumParts(tree)

    data <- rbind(
      data.table(
        node = tree$edge[, 2],
        regime = as.factor(PCMTreeGetRegimesForEdges(tree))),
      data.table(
        node = N+1L,
        regime = tree$part.regime[PCMTreeGetPartsForNodes(tree, N+1L)]))

    plotTree <- ggtree::`%<+%`(ggtree::ggtree(tree, ...), data)

    plotTree + aes(color = regime) +
      scale_color_manual(name = "regime", values = palette)
  } else {
    stop("ERR:026i0:PCMBase:PCMTree.R:PCMTreePlot:: Calling PCMTreePlot needs ggtree package to be installed from Bioconductor. Check the instructions at https://bioconductor.org/packages/release/bioc/html/ggtree.html. Ggtree was not on CRAN at the time of releasing PCMBase and is not declared as dependency in the PCMBase-description.")
  }
}
