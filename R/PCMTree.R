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
    stop("ERR:020d0:PCMBase:PCM.R:PCMTreeSetRegimes:: argument tree should be a phylo.")
  }
  if(!is.null(regimes)) {
    if(length(regimes) != length(nodes) + 1 ||
       length(unique(regimes)) < length(regimes)) {
      stop("ERR:020d1:PCMBase:PCM.R:PCMTreeSetRegimes:: regimes should be a character or integer vector of length equal to length(nodes) + 1 and not have duplicated elements.")
    }
  }
  preorder <- PCMTreePreorder(tree)
  edge.regime <- rep(1, nrow(tree$edge))
  nextRegime <- 2

  endingAt <- order(rbind(tree$edge, c(0, N+1))[, 2])

  for(ei in preorder) {
    i <- tree$edge[ei, 2]
    j <- tree$edge[ei, 1]
    ej <- endingAt[j]

    if(i %in% nodes) {
      edge.regime[ei] <- nextRegime
      nextRegime <- nextRegime + 1
    } else {
      edge.regime[ei] <- edge.regime[ej]
    }
  }
  edge.regime
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
#' @param regimes a vector of the same mode as tree$edge.regime (i.e. character or
#' integer) giving the order of the regimes on the tree in which the starting nodes
#' will be given. Default \code{PCMTreeUniqueRegimes(tree)}.
#' @return an integer with elements equal to the starting nodes for each regime in
#' \code{regimes}.
#' @seealso \code{\link{PCMTreeSetRegimes}}
#' @export
PCMTreeGetStartingNodesRegimes <- function(tree, regimes = PCMTreeUniqueRegimes(tree)) {
  if(!inherits(tree, "phylo")) {
    stop("ERR:020g0:PCMBase:PCM.R:PCMTreeGetStartNodesRegimes:: argument tree should be a phylo.")
  }

  if(!all(regimes %in% PCMTreeUniqueRegimes(tree))) {
    stop("ERR:020g1:PCMBase:PCM.R:PCMTreeGetStartNodesRegimes:: some of the regimes given are not
         found as regimes in the tree.")
  }

  N <- PCMTreeNumTips(tree)
  nodes <- integer(length(regimes))
  names(nodes) <- as.character(regimes)

  # start with all regimes being descending from the root; we will update them
  # as we encounter them in a preorder traversal.
  nodes[] <- 0

  preorder <- PCMTreePreorder(tree)

  for(ei in preorder) {
    if(as.character(tree$edge.regime[ei]) %in% regimes &&
       nodes[as.character(tree$edge.regime[ei])] == 0) {
      i <- tree$edge[ei, 2]
      j <- tree$edge[ei, 1]
      if(j != N+1) {
        nodes[as.character(tree$edge.regime[ei])] <- i
      } else {
        nodes[as.character(tree$edge.regime[ei])] <- N+1
      }
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
    stop("ERR:02010:PCMBase:PCM.R:PCMRegimesUniqueTree:: tree$edge.regime is NULL,
         but should be a character or an integer vector denoting regime names.")
  }
  uniqueRegimesPos <- integer(PCMTreeNumUniqueRegimes(tree))
  uniqueRegimes <- sort(unique(tree$edge.regime))
  names(uniqueRegimesPos) <- as.character(uniqueRegimes)
  uniqueRegimesPos[] <- 0L
  preorder <- PCMTreePreorder(tree)
  currentPos <- 1
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
PCMTreeJumps <- function(tree) {
  if(!is.null(tree$edge.jump)) {
    if(length(tree$edge.jump) != nrow(tree$edge) |
       !isTRUE(all(tree$edge.jump %in% c(0L, 1L)))) {
      stop("ERR:02081:PCMBase:PCM.R:PCMTreeJumps:: tree$edge.jump should
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
  regimes <- match(tree$edge.regime, PCMRegimes(model))
  if(any(is.na(regimes))) {
    stop(paste0("ERR:02071:PCMBase:PCM.R:PCMTreeMatchRegimesWithModel:: ",
                " Some of the regimes in tree$edge.regime not found in",
                "attr(model, 'regimes').\n",
                "unique regimes on the tree:", toString(PCMTreeUniqueRegimes(tree)), "\n",
                "attr(model, 'regimes')", toString(PCMRegimes(model))))
  }
  regimes
}


