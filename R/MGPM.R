# Copyright 2016-2020 Venelin Mitov
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

#' Construct a MGPM object from a vector
#'
#' MGPM stays for Mixed Gaussian Phylogenetic Model. A MGPM object represents a
#' list of all variables that are subject to change during a probabilistic
#' inference of such a model for a given tree and trait data (see Details). This
#' is unlike a \code{\link{MGPMContext}} object that contains the necessary
#' data and meta-information that remains constant (e.g. model types) or changes
#' rarely (e.g. metaI objects) during the inference.
#'
#' @param s a real (double) vector representing a MGPM. Some of the
#'  entries in this vector are coerced to integers denoting nodes in the
#'  tree, regimes or model types (see example and \code{\link{MGPMVector}}).
#' @param ctx a MGPMContext object.
#'
#' @details
#' A MGPM can be encoded as a numerical vector:
#' \deqn{\vec{s}=(K, n_2,...,n_{K+1}, l_2,...,l_{K+1}, r_2,...,r_{K+1}, m_1,...,m_R, v_1,...,v_P)^T}
#' Use the constructor \code{\link{MGPMVector}} to create such vectors. The
#' vector elements are described as follows:
#' \describe{
#' \item{\eqn{K}: }{number of shifts;}
#' \item{\eqn{(n_2,...,n_{K+1})^T}: }{shift nodes - the shifts in the model occur
#' at points within the branches leading to the shift nodes in tip-ward
#' direction. The corresponding locations of these points are specified by
#' \eqn{(l_2,...,l_{K+1})^T}. The nodes \eqn{(n_2,...,n_{K+1})^T} should be ordered
#' according to \code{\link{PCMTreePreorder}(ctx$tree)}.}
#' \item{\eqn{(l_2,...,l_{K+1})^T}: }{offset of the shift points measured as
#' distances from the beginnings of the branches leading to shift nodes in
#' tip-ward direction.}
#' \item{\eqn{(r_2,...,r_{K+1})^T}: }{regime index vector. This is an integer
#' vector with elements among \eqn{(N+1,n_2,...,n_{K+1})^T}, indicating the regime
#' associated with each part in the tree. The regimes are named as the shift
#' nodes, with N+1 corresponding to the part (and regime) originating at the
#' root. It is possible to have lumped regimes, that is, different parts
#' of the tree having the same regime. This regime-lumping must obey the
#' following rules:
#' \enumerate{
#' \item neighbor parts cannot have a lumped regime. Two parts originating at
#' nodes \eqn{n_i} and \eqn{n_j} in the tree are called neighbor parts if they
#' are separated solely by \eqn{n_i} or by \eqn{n_j};
#' \item to resolve the conflict between the shift nodes of the different parts
#' covered by a lumped regime, it is established that the name of a lumped
#' regime must equal the shift-node that appears first, according to
#' \code{\link{PCMTreePreorder}(ctx$tree)}.
#' }}
#' \item{\eqn{(m_1,...,m_R)^T}: }{model type assignment to the unique regimes.
#' This is an integer vector with elements between 1 and M, M denoting the
#' total number of model types possible. Each element corresponds to an
#' element in \code{unique(c(N+1,r_2,...,r_{K+1}))}}
#' \item{\eqn{(v_1,...,v_P)^T}: }{real numbers passed to
#' \code{\link{PCMParamLoadOrStore}}. This is a vectorized form of the model
#' parameters.}
#' }
#'
#' @return an object of S3 class 'MGPM'.
#' @examples
#' # The PCMBase package comes with a collection of simulated objects, which we
#' # can use as example.
#'
#' tree <- PCMBaseTestObjects$tree.ab
#' model <-
#'   PCMExtractDimensions(PCMBaseTestObjects$model_MixedGaussian_ab, dims = 1:2)
#' X <- PCMBaseTestObjects$traits.ab.123[1:2, ]
#'
#' # Create a MGPM context
#' ctx <- MGPMContext(
#'   X, tree,
#'   model = MixedGaussian(
#'     k = 2,
#'     modelTypes = MGPMDefaultModelTypes(),
#'     mapping = structure(1:6, names = LETTERS[1:6]),
#'     Sigmae_x = Args_MixedGaussian_MGPMDefaultModelTypes()$Sigmae_x))
#'
#' state <- MGPM(s = c(
#'   3,              # K: number of shifts
#'   52, 46, 73,     # n_2, ..., n_{K+1}: shift nodes
#'   0.1, 0.4, 0.3,  # l_2, ..., l_{K+1}: offsets of the shift points relative
#'                   # to the beginnings of shift branches in tip-ward direction.
#'   52, 46, 41,     # r_2,...,r_{K+1}: regime indices corresponding to the shifts.
#'                   # The regime of the part starting at the root node (41) is
#'                   # set to 41 (not included in the list). The regime for the
#'                   # part 73 is again 41, meaning that this regime is lumped with
#'                   # the root regime. So the number of regimes is R=3.
#'   5, 2, 2         # model type mapping for the three regimes.
#'   ), ctx)
#'
#' @seealso \code{\link{MGPMVector}} \code{\link{MGPMContext}}
#' \code{\link{MGPMPosK}} \code{\link{MixedGaussian}}
#' @export
MGPM <- function(s, ctx) {
  # number of tips
  N <- PCMTreeNumTips(ctx$tree)

  # Initialize default values:
  K <- 0L
  n <- integer(0L)
  l <- double(0L)
  r <- N+1L
  m <- 1L
  v <- NULL

  R <- 1L
  P <- NA_integer_

  model <- NULL
  # The node labels in ctx$tree must be 1,2,...,M, with N+1 being the root and
  # 1,...,N being the tips.
  tree <- ctx$tree

  if(length(s) >= 1L) {
    # number of shifts
    K <- as.integer(s[1L])
    if(K > 0L && !(length(s) >= 1L + K + K + K)) {
      stop(paste0("MGPM:: The vector s is wrong length, given that K=", K,
                  ". It should have at least ", 1L + K + K + K, " elements."))
    } else {
      # shift nodes
      n <- as.integer(s[1L + seq_len(K)])

      # shift locations on the branches leading to the shift nodes
      l <- s[1L + K + seq_len(K)]

      # regimes: these are integers among n
      r <- as.integer(s[1L + K + K + seq_len(K)])

      # order n, l and r according to PCMTreePreorder
      ordern <- order(match(n, ctx$tree$edge[ctx$preorderNodes, 2]))

      n <- n[ordern]
      l <- l[ordern]
      r <- r[ordern]

      regimeNames <- unique(c(N+1L, r))

      # Number of different regimes
      R <- length(regimeNames)

      minSingletonBranchLength <-
        getOption("PCMBayes.minSingletonBranchLength", 0.1)

      # Insert singleton nodes wherever needed
      needASingleton <- (l >= minSingletonBranchLength)

      if(sum(needASingleton) > 0L) {
        # ids of rows in tree$edge and tree$edge.length
        branchesNeedingSingletons <- match(n[needASingleton], tree$edge[, 2L])
        # Offsets of the singletons in root-ward direction
        posSingletons <-
          tree$edge.length[branchesNeedingSingletons] - l[needASingleton]

        if(sum(posSingletons > minSingletonBranchLength) > 0L)
          tree <- PCMTreeInsertSingletons(
            tree,
            n[needASingleton][posSingletons >= minSingletonBranchLength],
            posSingletons[posSingletons >= minSingletonBranchLength])
      }

      # After possible insertion of singleton nodes, the node indices are not
      # anymore the same, so we use the node labels to set the partition
      part.names <- as.character(c(N+1L, n))
      part.regime <- structure(character(length(part.names)), names = part.names)
      part.regime[] <- as.character(c(N+1L, r))
      PCMTreeSetPartRegimes(tree, part.regime, setPartition = TRUE)

      m <- rep(1L, R)
      # model types associated with the regimes
      if(length(s) >= 1L + K + K + K + R) {
        m <- as.integer(s[1L + K + K + K + seq_len(R)])
      }
      names(m) <- regimeNames

      namesModelTypes <-
        names(ctx$modelTemplate)[sapply(ctx$modelTemplate, is.PCM)]

      # Create the model (MixedGaussian)
      model <- do.call(
        MixedGaussian,
        c(list(k = ctx$k,
               modelTypes = ctx$modelTemplate[namesModelTypes],
               mapping = m),
          attr(ctx$modelTemplate, "spec")[
            setdiff(names(attr(ctx$modelTemplate, "spec")), namesModelTypes)]))

      # Load the parameter values into the model (the number of parameters is
      # stored in P)
      P <- as.integer(attr(model, "p"))
      v <- double(P)
      if(length(s) >= 1L + K + K + K + R + P) {
        v <- s[-seq_len(1L + K + K + K + R)]
      }
      P1 <- PCMParamLoadOrStore(model, v, offset = 0, k = ctx$k, R = R,load = TRUE)
      if(P != P1) {
        stop(paste0(
          "MGPM:: something is wrong with the attribute p of the model",
          " and the number of parameters returned by PCMParamLoadOrStore."))
      }
    }
  }

  structure(list(
    K = K, R = R, P = P, n = n, l = l, r = r, m = m,
    v = v, model = model, tree = tree),
    class = "MGPM")
}


#' Check if an object is of S3 class 'MGPM'
#' @param o an object.
#' @return a logical.
#' @export
is.MGPM <- function(o) {
  inherits(o, "MGPM")
}

#' If necessary, convert an object to an MGPM object.
#' @param s a \code{\link{MGPM}} or a \code{\link{MGPMVector}}
#' object.
#' @param ctx a \code{\link{MGPMContext}} object. This is needed only if
#'  \code{s} is \code{\link{MGPMVector}} object.
#'
#' @return an MGPM object corresponding to \code{s}. If \code{s} is already
#' a \code{\link{MGPM}} object, it is returned as is. Otherwise, a
#' \code{\link{MGPM}} object is constructed from \code{s} and \code{ctx}.
#' @export
as.MGPM <- function(s, ctx = NULL) {
  if(is.MGPM(s)) {
    s
  } else if(is.MGPMVector(s)) {
    if(!is.MGPMContext(ctx)) {
      stop(paste0(
        "as.MGPM: ctx must be a MGPMContext object if s is a ",
        "MGPMVector."))
    }
    s <- MGPM(s, ctx)
  } else {
    stop("as.MGPM: s be a MGPM or a MGPMVector object.")
  }
}

#' Construct a MGPM vector.
#'
#' @param K an integer denoting the number of shifts.
#' @param n an integer vector of length \code{K} denoting shift nodes.
#' @param l a real vector of length \code{K} denoting shift-point offsets in
#' tip-ward direction from the beginning of the branches leading to shift-nodes.
#' @param r an integer vector of length \code{K} denoting the regime
#' associated with each part (see Details in \code{\link{MGPM}}).
#' @param m an integer vector of length \code{R}, where \code{R} denotes the
#' total number of unique regimes (see Details in \code{\link{MGPM}}).
#' @param v a real vector of length \code{P}, where \code{P} denotes the total
#' number of numeric parameters of the MixedGaussian model (see Details in
#' \code{\link{MGPM}}).
#' @return a real (double) vector representing the concatenation of
#' \code{K, n, l, r, m, v}, with S3 class set to 'MGPMVector'.
#'
#' @seealso \code{\link{MGPM}}
#' @export
MGPMVector <- function(K, n, l, r, m, v) {
  structure(
    c(as.double(K), as.double(n), as.double(l), as.double(r), as.double(m),
      as.double(v)),
    class = c("MGPMVector"))
}

#' Check if an object is of S3 class 'MGPMVector'
#' @param o an object.
#' @return a logical.
#' @export
is.MGPMVector <- function(o) {
  inherits(o, "MGPMVector")
}


#' Context of a MGPM
#' @inheritParams PCMLik
#' @param model a template used to build MixedGaussian model objects.
#' @param K,n,l,r,m,v default (or initial) variables for the number of
#' shifts, the shift nodes, the shift-offsets, the regimes, the mapped model
#' types and the parameter values (see \code{\link{MGPM}}).
#' By default, these are
#' \code{K = 0L, n = integer(0), l = double(0), r = N+1L, m = 1L, v = double(P)}.
#' Here, N is the number of tips in the tree, and P is the number of numerical
#' parameters in a MixedGaussian with a single regime mapped to model type 1.
#'
#' @return an object of S3 class 'MGPMContext'.
#' @export
MGPMContext <- function(
  X, tree, model, SE = matrix(0, PCMNumTraits(model), PCMTreeNumTips(tree)),
  K = 0L, n = integer(0L), l = double(0L), r = N+1L, m = 1L, v = NULL) {

  ctx <- list(
    k = nrow(X),
    X = X,
    SE = SE,
    modelTemplate = model,
    treeOriginal = tree,
    tree = PCMTree(tree))
  PCMTreeSetLabels(ctx$tree)

  class(ctx) <- "MGPMContext"

  ctx$MGPMTemplate <- MGPM(MGPMVector(K, n, l, r, m, v), ctx)

  ctx
}

#' @title Indices of different parts of an MGPM in an MGPMVector
#'
#' @description The different parts of an MGPM are described in
#' \code{\link{MGPMVector}}. \code{posK} returns the osition of the
#'  number of shifts (this is always equal to 1). See the section functions
#'  for the others.
#'
#' @param s a \code{\link{MGPM}} or a \code{\link{MGPMVector}}
#' object.
#' @param ctx a \code{\link{MGPMContext}} object. This is needed only if
#'  \code{s} is \code{\link{MGPMVector}} object.
#'
#' @return an integer vector.
#'
#' @seealso \code{\link{MGPM}}, \code{\link{MGPMVector}}
#' @export
MGPMPosK <- function(s, ctx = NULL) {
  1L
}

#' @describeIn MGPMPosK
#'
#' Positions of the shift nodes;
#'
#' @export
MGPMPosn <- function(s, ctx = NULL) {
  s <- as.MGPM(s, ctx)
  1L + seq_len(s$K)
}

#' @describeIn MGPMPosK
#'
#' Positions of the offsets of the shift-points in tip-ward direction
#' relative to the beginnings of the branches leading to shift nodes;
#'
#' @export
MGPMPosl <- function(s, ctx = NULL) {
  s <- as.MGPM(s, ctx)
  1L + s$K + seq_len(s$K)
}

#' @describeIn MGPMPosK
#'
#' Position of the regime id's for each shift node;
#'
#' @export
MGPMPosr <- function(s, ctx = NULL) {
  s <- as.MGPM(s, ctx)
  1L + s$K + s$K + seq_len(s$K)
}

#' @describeIn MGPMPosK
#'
#' Positions of the model type id's for each regime;
#'
#' @export
MGPMPosm <- function(s, ctx = NULL) {
  s <- as.MGPM(s, ctx)
  1L + s$K + s$K + s$K + seq_len(s$R)
}

#' @describeIn MGPMPosK
#'
#' Positions of the model parameters.
#'
#' @export
MGPMPosv <- function(s, ctx = NULL) {
  s <- as.MGPM(s, ctx)
  1L + s$K + s$K + s$K + s$R + seq_len(s$P)
}
