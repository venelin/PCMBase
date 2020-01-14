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

#' Construct a MGPM object
#'
#' MGPM stays for Mixed Gaussian Phylogenetic Model. A MGPM object represents a
#' list of all variables that are subject to change during a probabilistic
#' inference of such a model for a given tree and trait data (see Details). This
#' is unlike a \code{\link{MGPMContext}} environment that contains the necessary
#' data and meta-information that remains constant (e.g. model types) or changes
#' rarely (e.g. metaI objects) during the inference.
#'
#' @param ctx a \code{\link{MGPMContext}} object.
#' @param K a single non-negative integer number denoting the number of shifts
#' in the model. Alternatively, this can be a double vector equal to the
#' concatenation of the double-converted arguments \code{K, n, l, r, m, v}. In
#' this case, the arguments \code{n, l, r, m, v} will be ignored. To check that
#' this is the case, the length of the argument K is checked for being bigger
#' than 1. This makes sense because the shortest possible concatenation of
#' \code{K, n, l, r, m, v} is of length 2. Default: \code{0L}. See also
#' Details for details on the definition of \code{K, n, l, r, m, v}.
#' @param n an integer vector of length \code{K} containing ids of nodes in
#' \code{ctx$tree}. Ignored if \code{length(K) > 1}. Default value:
#' \code{integer(0)}. See also Details.
#' @param l a non-negative double vector of length \code{K}. Ignored if
#' \code{length(K) > 1}. Default value: \code{integer(0)}. See also Details.
#' @param r an integer vector of length \code{K} containing regime ids. Ignored
#' if \code{length(K) > 1}. Default value: \code{integer(0)}. See also Details.
#' @param m an integer vector of length R, where R is the number of unique
#' regimes, specifying model type mapping for each regime. Ignored if
#' \code{length(K) > 1}. Default: \code{1L}. See also Details.
#' @param v a double vector of length P where P is the number of variable
#' parameters of the MixedGaussian model corresponding to this MGPM object.
#' Ignored if \code{length(K) > 1}. By default this is set to
#' \code{double(P))}. See also Details.
#'
#' @details
#' A MGPM can be encoded as a numerical vector:
#' \deqn{\vec{s}=(K, n_2,...,n_{K+1}, l_2,...,l_{K+1}, r_2,...,r_{K+1}, m_1,...,m_R, v_1,...,v_P)^T}
#' Use the constructor \code{\link{MGPMVector}} to create such vectors. The
#' vector elements are described as follows:
#' \describe{
#' \item{\eqn{K}: }{number of shifts;}
#' \item{\eqn{(n_2,...,n_{K+1})^T}: }{shift nodes - the shifts in the model
#' occur at points within the branches leading to the shift nodes in tip-ward
#' direction. The corresponding locations of these points are specified by
#' \eqn{(l_2,...,l_{K+1})^T}. The nodes \eqn{(n_2,...,n_{K+1})^T} should be
#' ordered in increasing order.}
#' \item{\eqn{(l_2,...,l_{K+1})^T}: }{offsets of the shift points measured as
#' distances from the beginnings of the branches leading to shift nodes in
#' tip-ward direction. }
#' \item{\eqn{(r_2,...,r_{K+1})^T}: }{regime index vector. This is an integer
#' vector with elements among \eqn{(N+1,n_2,...,n_{K+1})^T}, indicating the
#' regime associated with each part in the tree. The regimes are named as the
#' shift nodes, with N+1 corresponding to the part (and regime) originating at
#' the root. The regime \eqn{r_1 = N+1} is always present and, therefore,
#' omitted. It is possible to have lumped regimes, that is, different parts
#' of the tree having the same regime. This regime-lumping must obey the
#' following rules:
#' \enumerate{
#' \item neighbor parts cannot have a lumped regime. Two parts originating at
#' nodes \eqn{n_i} and \eqn{n_j} in the tree are called neighbor parts if they
#' are separated solely by \eqn{n_i} or by \eqn{n_j};
#' \item to resolve the conflict between the shift nodes of the different parts
#' covered by a lumped regime, it is established that the name of a lumped
#' regime must equal the smallest of the lumped shift-nodes}.
#' }
#' \item{\eqn{(m_1,...,m_R)^T}: }{model type assignment to the unique regimes.
#' This is an integer vector with elements between 1 and M, M denoting the
#' total number of model types possible. Each element corresponds to an
#' element in \code{sort(unique(c(N+1,r_2,...,r_{K+1})))}}
#' \item{\eqn{(v_1,...,v_P)^T}: }{real numbers passed to
#' \code{\link{PCMParamLoadOrStore}}. This is a vectorized form of the model
#' parameters.}
#' }
#'
#' @return an object of S3 class 'MGPM'.
#' @examples
#' # The PCMBase package comes with a collection of simulated objects, which we
#' # can use for the example.
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
#' mgpmDefault <- MGPM(ctx)
#'
#' mgpm <- MGPM(ctx,
#'   K = 4,                    # number of shifts
#'   n = c(2, 46, 52, 73),     # n_2, ..., n_{K+1}: shift nodes
#'   l = c(0, 0.4, 0.1, 0.05), # l_2, ..., l_{K+1}: offsets of the shift points
#'                             # relative to the beginnings of shift branches in
#'                             # tip-ward direction.
#'   r = c(2, 46, 52, 41),     # r_2,...,r_{K+1}: regime indices corresponding
#'                             # to the shifts. The regime of the part starting
#'                             # at the root node (41) is set to 41 (not
#'                             # included in the list). The regime for the part
#'                             # 73 is again 41, meaning that this regime is
#'                             # lumped with the root regime. So the number of
#'                             # regimes is R=4, although there are 5 different
#'                             # parts in the tree.
#'   m = c(1, 5, 2, 3)         # model type mapping for the R regimes.
#'                             # Note that these correspond to
#'                             # sort(unique(c(N+1L, r_2, ..., r_{K+1}))).
#'                             # Hence, model type 1 corresponds to the
#'                             # single-branch regime (ending at tip 2), model
#'                             # type 5 correspods to the root regime (41),
#'                             # model type 2 corresponds to 46 and model type 3
#'                             # corresponds to 52.
#'   )
#'
#'
#' mgpm$n
#'
#' # the names correspond to the shift nodes.
#' mgpm$l
#' mgpm$r
#'
#' # The names of mgpm$m denote the regimes, while the values denote the model
#' # types:
#' mgpm$m
#'
#' mgpmVec <- as.MGPMVector(mgpm)
#' mgpmVec[MGPMPosv(mgpm)] <- PCMParamRandomVecParams(mgpm$model)
#'
#' mgpm2 <- as.MGPM(mgpmVec, mgpm$ctx)
#'
#' stopifnot(identical(mgpm$n, mgpm2$n))
#' stopifnot(identical(mgpm$l, mgpm2$l))
#' stopifnot(identical(mgpm$r, mgpm2$r))
#' stopifnot(identical(mgpm$m, mgpm2$m))
#'
#' # These two should be different values but same length:
#' mgpm$v
#' mgpm2$v
#'
#' mgpm3 <- MGPM(ctx,
#'   c(4, 2, 46, 52, 73, 0, 0.4, 0.1, 0.05, 2, 46, 52, 41, 1, 5, 2, 3))
#' stopifnot(identical(mgpm$n, mgpm3$n))
#' stopifnot(identical(mgpm$l, mgpm3$l))
#' stopifnot(identical(mgpm$r, mgpm3$r))
#' stopifnot(identical(mgpm$m, mgpm3$m))
#'
#' @seealso \code{\link{MGPMVector}} \code{\link{MGPMContext}}
#' \code{\link{MGPMPosK}} \code{\link{MixedGaussian}}
#' @export
MGPM <- function(
  ctx,
  K = 0L, n = integer(0L), l = double(0L), r = integer(0L), m = 1L, v = NULL) {

  if(length(K) > 1L) {
    s <- as.double(K)
  } else {
    s <- NULL
  }

  # number of tips in the tree
  N <- PCMTreeNumTips(ctx$tree)

  R <- 1L
  P <- NA_integer_

  model <- NULL

  # The node labels in ctx$tree must be 1,2,...,M, with N+1 being the root and
  # 1,...,N being the tips.
  tree <- ctx$tree

  if(is.null(s)) {
    # do nothing
  } else {
    # the first element is the number of shifts
    K <- as.integer(s[1L])
  }


  if(is.null(s)) {
    if(!(length(n) == K && length(l) == K && length(r) == K)) {
      stop(
        paste0("MGPM:: Some of n, l, r are of length different than K(", K, ")."))
    }
    nlr <- data.table(
      n = c(N + 1L, as.integer(n)),
      l = as.double(c(0.0, l)),
      r = c(N + 1L, as.integer(r)))
  } else {
    nlr <- data.table(
      n = c(N + 1L, as.integer(s[1L + seq_len(K)])),
      l = as.double(c(0.0, s[1L + K + seq_len(K)])),
      r = c(N + 1L, as.integer(s[1L + K + K + seq_len(K)])))
  }

  # avoid "no visible binding NOTE":
  i <- .I <- r2 <- NULL
  setkey(nlr, n)
  nlr[, i:=.I]

  # Row index in nlr corresponding to the root.
  iRoot <- nlr[n == N+1L, i]

  # n, l and r should be ordered in the increasing order of n
  ordern <- nlr[-iRoot, order(n)]
  if(!identical(ordern, seq_len(K))) {
    stop("MGPM:: the shift nodes should be ordered in increasing order.")
  }


  # Check that lumped regimes and set their names to the smallest shift-node
  nlr[, r2:=min(n), by=r]
  if(!nlr[, identical(r, r2)]) {
    print(nlr)
    stop(paste0(
      "MGPM:: each regime lumping several parts of the tree should be named ",
      "as the smallest shift-node among the shift-nodes for the lumped parts."))
  }
  if(length(nlr[, unique(c(n, r))]) != K+1L) {
    print(nlr)
    stop(paste0(
      "MGPM:: There are duplicated root or shift-nodes or some of the ",
      "regime names may not be among the root and the shift-nodes. \n",
      "root node and shift nodes: ", toString(nlr[, n]), "\n",
      "regimes: ", toString(nlr[, r]), "."))
  }

  # unique regime names
  regimeNames <- nlr[, unique(r)]

  # Number of unique regimes
  R <- length(regimeNames)

  minSingletonBranchLength <-
    getOption("PCMBayes.minSingletonBranchLength", 0.1)

  # Insert singleton nodes wherever needed
  needASingleton <- (nlr[, l] >= minSingletonBranchLength)

  if(sum(needASingleton) > 0L) {
    # ids of rows in tree$edge and tree$edge.length
    branchesNeedingSingletons <- nlr[, match(n[needASingleton], tree$edge[, 2L])]
    # Offsets of the singletons in root-ward direction
    posSingletons <-
      tree$edge.length[branchesNeedingSingletons] - nlr[, l[needASingleton]]

    if(sum(posSingletons > minSingletonBranchLength) > 0L)
      tree <- PCMTreeInsertSingletons(
        tree,
        nlr[, n[needASingleton][posSingletons >= minSingletonBranchLength]],
        posSingletons[posSingletons >= minSingletonBranchLength])
  }

  # After possible insertion of singleton nodes, the node indices are not
  # anymore the same, so we use the node labels to set the partition
  part.names <- nlr[, as.character(n)]
  part.regime <- structure(nlr[, as.character(r)], names = part.names)
  PCMTreeSetPartRegimes(tree, part.regime, setPartition = TRUE)

  if(is.null(s)) {
    if(length(m) != R) {
      stop(paste0(
        "MGPM:: the length of m (", length(m),
        ") should equal the number of unique regimes (", R, ")"))
    }
    m <- as.integer(m)
  } else {
    # model types associated with the regimes
    if(length(s) >= 1L + K + K + K + R) {
      m <- as.integer(s[1L + K + K + K + seq_len(R)])
    } else {
      stop(
        paste0("MGPM:: not enough model types supplied in input (s); using",
               "default model mapping 1 for all regimes."))
    }
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

  if(is.null(s)) {
    if(!is.null(v) && length(v) != P) {
      stop(paste0(
        "MGPM:: if v is specified, the length of v (", length(v),
        ") should equal the number of variable parameters of the model (", P, ")"))
    } else {
      v <- double(P)
    }
  } else {
    v <- double(P)
    if(length(s) >= 1L + K + K + K + R + P) {
      v <- s[-seq_len(1L + K + K + K + R)]
    }
    P1 <- PCMParamLoadOrStore(model, v, offset = 0, k = ctx$k, R = R,load = TRUE)
    if(P != P1) {
      stop(paste0(
        "MGPM:: something is wrong with the attribute p of the model ",
        "and the number of parameters returned by PCMParamLoadOrStore. ",
        "This could be a bug."))
    }
  }

  structure(list(
    K = K, R = R, P = P,
    n = nlr[-iRoot, n],
    l = structure(nlr[-iRoot, l], names = nlr[-iRoot, n]),
    r = structure(nlr[-iRoot, r], names = nlr[-iRoot, n]),
    m = m,
    v = v, model = model, tree = tree, ctx = ctx),
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
#' This is an S3 generic function.
#' @param o an object.
#' @param ctx a \code{\link{MGPMContext}} object. Default \code{NULL}.
#'
#' @return an MGPM object corresponding to \code{o}. If \code{o} is already
#' a \code{\link{MGPM}} object, it is returned as is. Otherwise, a
#' \code{\link{MGPM}} object is constructed from \code{o} and,
#' optionally \code{ctx}.
#' @export
as.MGPM <- function(o, ctx = NULL) {
  UseMethod("as.MGPM", o)
}

#' @export
as.MGPM.default <- function(o, ctx = NULL) {
  stop(paste0("as.MGPM: no method defined for class ", toString(class(o)), "." ))
}

#' @export
as.MGPM.MGPM <- function(o, ctx = NULL) {
  o
}

#' @export
as.MGPM.MGPMVector <- function(o, ctx = NULL) {
  if(!is.MGPMContext(ctx)) {
    stop(paste0(
      "as.MGPM: ctx must be a MGPMContext object if o is a ",
      "MGPMVector."))
  }
  MGPM(ctx, o)
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


#' Convert a MGPM object to a MGPMVector object
#' @param o a MGPMVector or a MGPM object
#' @return a MGPMVector object corresponding to o.
#' @seealso \code{\link{MGPM}} \code{\link{MGPMVector}}
#' @export
as.MGPMVector <- function(o) {
  if(is.MGPMVector(o)) {
    o
  } else if(is.MGPM(o)) {
    MGPMVector(o$K, o$n, o$l, o$r, o$m, o$v)
  } else {
    stop("as.MGPMVector: o must be a MGPM object.")
  }
}


#' Context of a MGPM
#' @inheritParams PCMLik
#' @param model a MixedGaussian model used as a template to build MixedGaussian
#' model objects.
#' @param tipModelType NULL or a character string indicating a member of
#' \code{model}. This argument allows to set a specific regime for all terminal
#' branches, such as a white noise model.
#' @param ... any additional objects to be included in the MGPMContext environment.
#'
#' @return an object of S3 class 'MGPMContext'.
#' @export
MGPMContext <- function(
  X, tree, model, SE = matrix(0, PCMNumTraits(model), PCMTreeNumTips(tree)),
  tipModelType = NULL, ...) {

  ctx <- list2env(list(
    k = nrow(X),
    X = X,
    SE = SE,
    tipModelType = tipModelType,
    modelTemplate = model,
    treeOriginal = tree,
    tree = PCMTree(tree)),
    ...)

  PCMTreeSetLabels(ctx$tree)

  class(ctx) <- "MGPMContext"

  ctx$MGPMTemplate <- MGPM(
    ctx,
    MGPMVector(
      K = 0L, n = integer(0L), l = double(0L), r = integer(0L), m = 1L, v = NULL))

  ctx
}

#' Check if an object is of S3 class 'MGPMContext'
#' @param o an object.
#' @return a logical.
#' @export
is.MGPMContext <- function(o) {
  inherits(o, "MGPMContext")
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

