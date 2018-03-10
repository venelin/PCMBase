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

#' @name PCMBase
#'
#' @title A General Framework for Gaussian Phylogenetic Comparative Models
#'
#' @description
#'
#'
NULL

#' @name PCM
#' @title PCM S3 objects
#' @param model This argument can take one of the following forms:
#' \item a character vector of the S3-classes of the model object to be
#' created (one model object can have one or more S3-classes);
#' \item a list of character vectors as above denoting the S3 classes of each
#' regime in a generalized multiple regime model. If available, the names of the
#' list-members will be used as regime names;
#' \item an S3 object representing a PCM model
#' @return an object of S3 class as defined by the model argument.
#' @export
PCM <- function(model, k = 1, regimes = 1,
                params = NULL, vecParams = NULL, offset = 0, ...) {
  UseMethod("PCM", model)
}

#' @describeIn PCM
#' @export
PCM.default <- function(model, k = 1, regimes = 1,
                        params = NULL, vecParams = NULL, offset = 0, ...) {
  stop(paste0("ERR:02091:PCMBase:MultivariatePCM.R:PCM.default:: You should provide a PCM object, but provided a ", class(model)[1]))
}

#' Create a PCM
#' @inheritParams PCM
#' @export
PCM.PCM <- function(model, k = 1, regimes = 1,
                    params = NULL, vecParams = NULL, offset = 0, ...) {

  if(is.null(attr(model, "description", exact = TRUE))) {
    attr(model, "description") <- PCMDescribe(model, ...)
  }

  if(is.null(attr(model, "specParams", exact = TRUE))) {
    attr(model, "specParams") <- PCMSpecifyParams(model, ...)
  }

  if(is.null(attr(model, "p", exact = TRUE))) {
    attr(model, "p") <- PCMNumParams(model)
  }

  specParams <- attr(model, "specParams")

  # TODO: parameter validation checks
  for(name in names(specParams)) {
    model[[name]] <- specParams[[name]][["default"]]
  }

  if(!is.null(params)) {
    PCMSetParams(model, params, ...)
  }

  if(!is.null(vecParams)) {
    PCMSetOrGetVecParams(model, vecParams, offset, ...)
  }

  model
}

#' @describeIn PCM
#' @export
PCM.character <- function(model, k = 1, regimes = 1,
                          params = NULL, vecParams = NULL, offset = 0, ...) {
  modelObj <- list()
  if(! (model %in% "PCM") ) {
    model <- c(model, "PCM")
  }
  class(modelObj) <- model
  attr(modelObj, "k") <- k
  attr(modelObj, "regimes") <- regimes

  PCM(modelObj, k, regimes, params, vecParams, offset, ...)
}

#' @describeIn PCM
#' @export
is.PCM <- function(x) inherits(x, "PCM")

#' @describeIn PCM
#' @export
print.PCM <- function(x, ...) cat (format(x, ...), "\n")

#' @describeIn PCM
#' @export
format.PCM <- function(x, ...) {
  if(!is.PCM(x)) {
    stop("ERR:02061:PCMBase:MultivariatePCM.R:format.PCM:: x must inherit from S3 class PCM.")
  }
  res <- paste0("PCM of S3 class ", class(x)[1],
                "; k=", attr(x, "k"),
                "; p=", attr(x, "p"),
                "; regimes: ", toString(attr(x, "regimes")), ". Parameters:\n")
  for(name in names(x)) {
    strList <- as.list(capture.output(print(x[[name]])))
    strList$sep = "\n"
    res <- paste0(res, name, ":\n", do.call(paste, strList))
    res <- paste0(res, "\n")
  }
  res
}

#' Human friendly description of a PCM
#' @param model a PCM model object
#' @details This S3 generic function is intended to be specified for user models
#' @return a character string
#' @export
PCMDescribe <- function(model, ...) {
  UseMethod("PCMDescribe", model)
}

#' @describeIn PCMDescribe
#' @export
PCMDescribe.PCM <- function(model, ...) {
  "PCM base class model with no parameters; serves as a basis for PCM model classes"
}

#' Parameter specification of PCM model
#' @param model a PCM model object
#' @return a list
#' @export
PCMSpecifyParams <- function(model, ...) {
  UseMethod("PCMSpecifyParams", model)
}

#' @describeIn PCMSpecifyParams
#' @export
PCMSpecifyParams.PCM <- function(model, ...) {
  list()
}

#' Number of traits modeled by a PCM
#' @param model a PCM object
#' @return an integer
#' @export
PCMNumTraits <- function(model) {
  UseMethod("PCMNumTraits", model)
}

#' @describeIn PCMNumTraits
#' @export
PCMNumTraits.PCM <- function(model) {
  attr(model, "k")
}

#' Number of free parameters describing fully a PCM
#' @param model a PCM object
#' @param ... other arguments (possible future use)
#' @return an integer
#' @export
PCMNumParams <- function(model, ...) {
  UseMethod("PCMNumParams", model)
}

#' @describeIn PCMNumParams
#' @export
PCMNumParams.PCM <- function(model, ...) {
  k <- attr(model, "k")
  R <- length(attr(model, "regimes"))
  specParams <- attr(model, "specParams")

  vecParamSizes <- c(rep = 1, full = k)
  matParamSizes <- c(diag1 = 1, diag = k, upper.tri = .5*k*(k-1), upper.tri.diag = .5*k*(k+1),
                     lower.tri = .5*k*(k-1), lower.tri.diag = .5*k*(k+1), symmetric = .5*k*(k+1), full=k*k)
  p <- 0
  for(name in names(specParams)) {
    type <- specParams[[name]]$type
    indices <- specParams[[name]]$indices

    if(type[1] == "model") {
      p <- p + PCMNumParams(model[[name]], ...)
    } else if(type[1] == "gvector") {
      if(type[2] == "custom") {
        ind <- indices(p, k)
        p <- p + length(unique(ind[ind>p]))
      } else {
        p <- p + vecParamSizes[type[2]]
      }
    } else if(type[1] == "gmatrix") {
      if(type[2] == "custom") {
        ind <- indices(p, k)
        p <- p + length(unique(ind[ind>p]))
      } else {
        p <- p + matParamSizes[type[2]]
      }
    } else if(type[1] == "vector") {
      if(type[2] == "custom") {
        for(r in 1:R) {
          ind <- indices(p, k)
          p <- p + length(unique(ind[ind>p]))
        }
      } else {
        p <- p + R*vecParamSizes[type[2]]
      }
    } else if(type[1] == "matrix") {
      if(type[2] == "custom") {
        for(r in 1:R) {
          ind <- indices(p, k)
          p <- p + length(unique(ind[ind>p]))
        }
      } else {
        p <- p + R*matParamSizes[type[2]]
      }
    }
  }
  p
}

#' Wrapper for length(tree$tip.label)
#' @param tree a phylo object
#' @return the number of tips in tree
#' @export
PCMNumTips <- function(tree) {
  length(tree$tip.label)
}

#' Number of all nodes in a tree
#'
#' @details Wrapper for nrow(tree$edge) + 1
#' @param tree a phylo object
#' @return the number of nodes in tree including root, internal and tips.
#' @export
PCMNumNodes <- function(tree) {
  nrow(tree$edge) + 1
}

#' Number of regimes on a tree modeled with PCM
#' @param tree a phylo object
#' @return the number of different regimes encountered on the tree branches
#' @export
PCMNumRegimes <- function(tree) {
  # number unique regimes
  if(is.null(tree$edge.regime)) {
    if(!is.null(names(tree$edge.length))) {
      tree$edge.regime <- names(tree$edge.length)
      regimesUniqueTree <- unique(tree$edge.regime)
    } else {
      regimesUniqueTree <- 1
    }
  } else {
    regimesUniqueTree <- unique(tree$edge.regime)
  }

  length(regimesUniqueTree)
}

#' Jumps in modeled traits associated with branches in a tree
#' @inheritParams PCMNumTips
#' @return an integer vector of 0's and 1's with entries correspondin to the
#' denoting if a jump took place at the beginning of a branch.
PCMJumps <- function(tree) {
  if(!is.null(tree$edge.jump)) {
    if(length(tree$edge.jump) != nrow(tree$edge) |
       !isTRUE(all(tree$edge.jump %in% c(0L, 1L)))) {
      stop("ERR:02081:PCMBase:MultivariatePCM.R:PCMJumps:: tree$edge.jump should
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
#'   denoting the regime on each branch in the tree
#' @export
PCMRegimes <- function(tree, model) {
  R <- PCMNumRegimes(tree)
  if(R == 1) {
    regimes <- rep(1, length(tree$edge.length))
  } else {
    regimes <- match(tree$edge.regime, attr(model, "regimes"))
    if(any(is.na(regimes))) {
      stop(paste0("ERR:02071:PCMBase:MultivariatePCM.R:PCMRegimes:: ",
                  " Some of the regimes in tree$edge.regime not found in",
                  "attr(model, 'regimes').\n",
                  "tree$edge.regime=", toString(tree$edge.regime), "\n",
                  "attr(model, 'regimes')", toString(attr(model, 'regimes'))))
    }
  }
  regimes
}



#' Load/store a vector parameter from/to a vector of all parameters in a model.
#'
#' @details This function has both, a returned value and side effects. By default the function
#' loads elements from vecParams into v. This behavior is reversed if the argument load
#' is set to FALSE.
#'
#' @param v numeric k-vector. If load==TRUE (default) this has to be a vector-object in
#' the calling environment (not a temporary object such as the result from an algebraic expression).
#' @param vecParams numeric vector containing all parameters in the model. If load==FALSE, this
#' has to be a vector in the calling environment (not a temporary object such as the result from
#' an algebraic expression).
#' @param offset integer denoting the offset in vecParams (0 means start from vecParams[1] onwards).
#' @param k integer denoting the dimensionality (length) of the resulting vector.
#' @param type a character string indicating the type of loading. possible values are:
#' \item{"rep"}{repeat \code{vecParams[offset+1]} k times}
#' \item{"full"}{copy \code{vecParams[offset + (1:k)]}}
#' \item{"custom"}{use a custom template vector which elements specified by mask are
#' to be replaced with \code{vecParams[indices(offset, k)]}}.
#' @param mask a logical vecgtor pointing to elements in template to be replaced. used only with type=="custom".
#' @param indices a function of the form function(offset, k) returning an integer vector of
#' length equal to the number of TRUE elements in mask, indicating the position in vecParams to
#' be copyed from. used only with type=="custom".
#' @param load a logical indicating if loading from or storing to vecParams should be done.
#' @return an integer denoting the number of elemnents read or written from/to vecParams.
#' In the case of type=="custom",
#' the number of indices bigger than offset returned by the function indices(offset, k).
#' @export
PCMLoadVectorParameter <- function(
  v, vecParams, offset, k,
  type = c("rep", "full", "custom"),
  mask = rep(TRUE, k), indices = function(offset, k) offset + (1:k),
  load = TRUE) {

  if(load) {
    "%op%" <- `<-`
    maskRep <- 1:k
  } else {
    "%op%"<- function(a,b) eval(substitute(b<-a), parent.frame())
    maskRep <- 1
  }

  if(type[1] == "rep") {
    num <- 1
    eval(substitute(v[maskRep] %op% vecParams[offset + 1]), parent.frame())
  } else if(type[1] == "full") {
    num <- k
    eval(substitute(v[] %op% vecParams[offset + (1:k)]), parent.frame())
  } else if(type[1] == "custom") {
    if(is.function(indices)) {
      ind <- indices(offset, k)
      num <- length(unique(ind[ind>offset]))
      maskCustom <- mask
      eval(substitute(v[maskCustom] %op% vecParams[ind]), parent.frame())
    } else {
      stop("ERR:020a1:PCMBase:MultivariatePCM.R:PCMLoadVectorParameter:: indices should be a
           function(offset, k) returning an integer vector.")
    }
  } else if(type[1] == "fixed") {
    num <- 0
  } else {
    stop(paste0("ERR:020a2:PCMBase:MultivariatePCM.R:PCMLoadVectorParameter:: type ", type[1], " not recognized."))
  }
  num
}

#' Load/store a matrix parameter from/to a vector of all parameters in a model.
#'
#' @details This function has both, a returned value and side effects. By default the function
#' loads elements from vecParams into m. This behavior is reversed if the argument load
#' is set to FALSE.
#'
#' @param m numeric k x k matrix. If load==TRUE (default) this has to be a matrix-object in
#' the calling environment (not a temporary object such as the result from an algebraic expression).
#' @param vecParams numeric vector containing all parameters in the model. If load==FALSE, this
#' has to be a vector in the calling environment (not a temporary object such as the result from
#' an algebraic expression).
#' @param k integer denoting the dimensionality (number of rows and colums) of the resulting matrix.
#' @param offset integer denoting the offset in vecParams (0 means start from vecParams[1] onwards).
#' @param type a character string indicating the type of loading. possible values are:
#' \item{"diag1"}{create a diagonal matrix with all diagonal elements repeating \code{vecParams[offset+1]}.}
#' \item{"diag"}{create a diagonal matrix with a diagonal copied from \code{vecParams[offset + (1:k)]}.}
#' \item{"upper.tri"}{load an upper triangular matrix with zero diagonal.}
#' \item{"upper.tri.diag"}{load an upper triangular matrix including diagonal.}
#' \item{"lower.tri"}{load a lower triangular matrix with zero diagonal.}
#' \item{"lower.tri.diag"}{load a lower triangular matrix including diagonal.}
#' \item{"symmetric"}{load a symmetric matrix. Only the elements of the upper triangle and diagonal are
#' loaded from vecParams.}
#' \item{"full"}{load a full.}
#' \item{"custom"}{use a custom template matrix which's elements specified by mask are
#' to be replaced with \code{vecParams[offset + indices]}}.
#' @param mask a logical k x k matrix pointing to elements in template to be replaced. used only with type=="custom".
#' @param indices a function of the form function(offset, k) returning an integer vector of
#' length equal to the number of TRUE elements in mask, indicating the position in vecParams to
#' be copyed from. used only with type=="custom".
#' @param load a logical indicating if loading from or storing to vecParams should be done.
#'
#' @return an integer denoting the number of elemnents read or written from/to vecParams.
#' In the case of type=="custom",
#' the number of indices bigger than offset returned by the function indices(offset, k).
#' @export
PCMLoadMatrixParameter <- function(
  m, vecParams, offset, k,
  type = c("diag1", "diag", "upper.tri", "upper.tri.diag", "lower.tri", "lower.tri.diag", "symmetric", "full", "custom", "fixed"),
  mask = matrix(TRUE, k, k),
  indices = function(offset, k) offset + (1:(k*k)),
  load = TRUE) {

  if(load) {
    "%op%" <- `<-`
    maskDiag1 <- 1:k
  } else {
    "%op%" <- function(a,b) eval(substitute(b<-a), parent.frame())
    maskDiag1 <- 1
  }

  if(type[1] == "diag1") {
    num <- 1
    eval(substitute(diag(m)[maskDiag1] %op% vecParams[offset + (1:num)]), parent.frame())
  } else if(type[1] == "diag") {
    num <- k
    eval(substitute(diag(m) %op% vecParams[offset + (1:num)]), parent.frame())
  } else if(type[1] == "upper.tri") {
    num <- k*(k-1)/2
    eval(substitute(m[upper.tri(m)] %op% vecParams[offset + (1:num)]), parent.frame())
  } else if(type[1] == "upper.tri.diag") {
    num <- k*(k+1)/2
    eval(substitute(m[upper.tri(m, diag=TRUE)] %op% vecParams[offset + (1:num)]), parent.frame())
  } else if(type[1] == "lower.tri") {
    num <- k*(k-1)/2
    eval(substitute(m[lower.tri(m)] %op% vecParams[offset + (1:num)]), parent.frame())
  } else if(type[1] == "lower.tri.diag") {
    num <- k*(k+1)/2
    eval(substitute(m[lower.tri(m, diag=TRUE)] %op% vecParams[offset + (1:num)]), parent.frame())
  } else if(type[1] == "symmetric") {
    num <- k*(k+1)/2
    eval(substitute({
      m[upper.tri(m, diag=TRUE)] %op% vecParams[offset + (1:num)]
      }), parent.frame())
    if(load) {
      eval(substitute({
        m[lower.tri(m)] <- 0
        m <- m+t(m)
        diag(m) <- 0.5 * diag(m)
      }), parent.frame())
    }
  } else if(type[1] == "full") {
    num <- k*k
    eval(substitute(m[,] %op% vecParams[offset+(1:num)]), parent.frame())
  } else if(type[1] == "custom") {
    if(is.function(indices)) {
      ind <- indices(offset, k)
      num <- length(unique(ind[ind>offset]))
      maskCustom <- mask
      eval(substitute(m[maskCustom] %op% vecParams[ind]), parent.frame())
    } else {
      stop("ERR:020b1:PCMBase:MultivariatePCM.R:PCMLoadMatrixParameter:: indices should be a
           function(offset, k) returning an integer vector.")
    }
  } else if(type[1] == "fixed") {
    num <- 0
  } else {
    stop(paste0("ERR:020b2:PCMBase:MultivariatePCM.R:PCMLoadMatrixParameter:: type ", type[1], " not recognized."))
  }
  num
}

#' Set model parameters from a named list
#' @param tree a phylo object (possible future use)
#' @param model a PCM model object
#' @param params a named list with elements among the names found in attr(model, "specParams")
#' @param inplace logical indicating if the parameters should be set "inplace" for the
#' model object in the calling environment or a new model object with the parameters set
#' as specified should be returned. Defaults to TRUE
#' @return If inplace is TRUE, the function only has a side effect of setting the
#' parameters of the model object in the calling environment; otherwise the function
#' returns a modified copy of the model object.
#' @export
PCMSetParams <- function(model, params, inplace = TRUE, ...) {
  UseMethod("PCMSetParams", model)
}

#' @describeIn PCMSetParams
#' @export
PCMSetParams.PCM <- function(model, params, inplace = TRUE, ...) {
  specParams <- attr(model, "specParams")
  for(name in names(params)) {
    if(! (name%in%names(specParams)) ) {
      stop(paste0("ERR:020c1:PCMBase:MultivariatePCM.R:PCMSetParams:: ", name,
                  " is not a settable parameter of the model, check attr(model, 'specParams')."))
    }
    type <- specParams[[name]]$type

    if(type != "model") {
      if(! identical(length(model[[name]]), length(params[[name]])) ) {
        stop(paste0("ERR:020c2:PCMBase:MultivariatePCM.R:PCMSetParams:: params[[", name,
                    "]] is not the same length as model[[", name, "]]; ",
                    "length(params[[", name, "]])=", length(params[[name]]),
                    ", length(model[[", name, "]])=", length(model[[name]]), "."))
      }
      if(! identical(dim(model[[name]]), dim(params[[name]])) ) {
        stop(paste0("ERR:020c3:PCMBase:MultivariatePCM.R:PCMSetParams:: params[[", name,
                    "]] is not the same dimension as model[[", name, "]]; ",
                    "dim(params[[", name, "]])=", str(dim(params[[name]])),
                    ", length(model[[", name, "]])=", str(dim(model[[name]])), "."))
      }
    }
  }
  for(name in names(params)) {
    if(inplace) {
      eval(substitute(model[[name]] <- params[[name]]), parent.frame())
    } else {
      model[[name]] <- params[[name]]
    }
  }

  if(!inplace) {
    model
  }
}

#' Inplace set or get the parameters of a PCM from or into a numeric vector
#'
#' @export
PCMSetOrGetVecParams <- function(
  model, vecParams, offset = 0, set = TRUE, ...) {
  UseMethod("PCMSetOrGetVecParams", model)
}

#' Set or get the parameters of a BM model from or into a numeric vector
#' @inheritParams PCMSetOrGetVecParams
#'
#' @export
PCMSetOrGetVecParams.PCM <- function(
  model, vecParams, offset = 0, set = TRUE, ... ) {

  p <- 0

  specParams <- attr(model, "specParams")
  k <- attr(model, "k")

  for(name in names(specParams)) {
    type <- specParams[[name]]$type
    mask <- specParams[[name]]$mask
    indices <- specParams[[name]]$indices


    if(type[1] == "vector") {
      for(r in 1:length(attr(model, "regimes"))) {
        p <- p + eval(substitute(PCMLoadVectorParameter(
          model[[name]][, r], vecParams, offset + p, k, type = type[2], mask = mask, indices = indices, load = set)),
          parent.frame())
      }
    } else if(type[1] == "matrix") {
      for(r in 1:length(attr(model, "regimes"))) {
        p <- p + eval(substitute(PCMLoadMatrixParameter(
          model[[name]][,, r], vecParams, offset + p, k, type = type[2], mask = mask, indices = indices, load = set)),
          parent.frame())
      }
    } else if(type[1] == "gvector") {
      # global vector for all regimes
      p <- p + eval(substitute(PCMLoadVectorParameter(
        model[[name]], vecParams, offset + p, k, type=type[2], mask = mask, indices = indices, load = set)),
        parent.frame())
    } else if(type[1] == "gmatrix") {
      # global vector for all regimes
      p <- p + eval(substitute(PCMLoadMatrixParameter(
        model[[name]], vecParams, offset + p, k, type=type[2], mask = mask, indices = indices, load = set)),
        parent.frame())
    } else if(type[1] == "model") {
      # nested PCM corresponding to a regime
      p <- p + eval(substitute(PCMSetOrGetVecParams(
        model[[name]], vecParams, offset + p, set, ...)), parent.frame())
    }
  }
  p
}

#' Simulation of a phylogenetic comparative model on a tree
#'
#' @description Generate trait data on a tree according to a multivariate stochastic
#' model with one or several regimes
#'
#' @param X0 a numeric vector of length k (the number of traits) specifying the
#' trait values at the root of the tree.
#' @param tree a phylo object specifying a rooted tree.
#' @param model an S3 object specifying the model (see Details).
#' @param metaI a named list containg meta-information about the data and the
#' model.
#'
#' @details Internally, this function uses the \code{\link{PCMCond}} iimplementation
#'  for the given model class.
#'
#' @return a list with two members as follows:
#' \item{values}{numeric M x k matrix of values at all nodes of the tree (root,
#' internal and tip), where M is the number of nodes: M=dim(tree$edge)[1]+1,
#' with indices from 1 to N=length(tree$tip.label) corresponding to tips, N+1
#' corresponding to the root and bigger than N+1 corresponding to internal nodes.}
#' \item{errors}{numeric M x k matrix of errors at all nodes of the tree (root,
#' internal and tip). Note that these errors are only simulated if the model contains
#' a k x k x R matrix member with the name 'Sigmae'. In this case the errors are
#' simulated as a random Gaussian noise with 0 mean and variance covariance matrix
#' defined by Sigmae. Note also that the errors on the root- and the internal nodes
#' are not used as input for the values/errors at their descending nodes, because
#' they are treated as non-heritable components or measurement error.}
#' The function will fail in case that the length of the argument vector X0 differs
#' from the number of traits specified in \code{metaI$k}. Error message:
#' "ERR:02002:PCMBase:MultivariatePCM.R:PCMSim:: X0 must be of length ...".
#'
#' @importFrom mvtnorm rmvnorm
#' @seealso \code{\link{PCMLik}} \code{\link{PCMValidate}} \code{\link{PCMCond}}
#' @export
PCMSim <- function(
  tree, model, X0,
  metaI = PCMInfo(X = NULL, tree = tree, model = model, verbose = verbose),
  verbose = FALSE) {

  if(length(X0)!=metaI$k) {
    stop(paste('ERR:02002:PCMBase:MultivariatePCM.R:PCMSim:: X0 must be of length', metaI$k, '.'))
  }

  values <- errors <- matrix(0, nrow=metaI$k, ncol=dim(tree$edge)[1]+1)
  values[, metaI$N + 1] <- X0

  ordBF <- metaI$preorder

  # create a list of random generator functions for each regime
  PCMCondObjects <- lapply(1:metaI$R, function(r) {
    PCMCond(tree, model = model, r = r, metaI = metaI, verbose = verbose)
  })

  for(edgeIndex in ordBF) {
    obj <- PCMCondObjects[[metaI$r[edgeIndex]]]
    if(!is.null(obj$random)) {
      values[, tree$edge[edgeIndex, 2]] <-
        obj$random(n=1, x0 = values[, tree$edge[edgeIndex,1]], t = tree$edge.length[edgeIndex], edgeIndex = edgeIndex)
    } else {
      values[, tree$edge[edgeIndex, 2]] <-
        PCMCondRandom(obj, n=1, x0 = values[, tree$edge[edgeIndex,1]], t = tree$edge.length[edgeIndex], edgeIndex = edgeIndex)
    }

    if(!is.null(model$Sigmae)) {
      errors[, tree$edge[edgeIndex, 2]] <-
        rmvnorm(1, rep(0, metaI$k),
                as.matrix(model$Sigmae[,, metaI$r[edgeIndex]]))
    }
  }

  list(values=values, errors=errors)
}

#' Likelihood of a multivariate Gaussian phylogenetic comparative model with non-interacting lineages
#'
#' @description The likelihood of a PCM represets the probability density function
#'   of observed trait values (data) at the tips of a tree given the tree and
#'   the model parameters. Seen as a function of the model parameters, the
#'   likelihood is used to fit the model to the observed trait data and the
#'   phylogenetic tree (which is typically inferred from another sort of data, such
#'   as an alignment of genetic sequences for the species at the tips of the tree).
#'   The \code{\link{PCMLik}} function
#'   provides a common interface for calculating the (log-)likelihood of different
#'   PCMs.
#'   Below we denote by N the number of tips, by M the total number of nodes in the
#'   tree including tips, internal and root node, and by k - the number of traits.
#'
#' @param X a \code{k x N} numerical matrix with possible \code{NA} entries. Each
#'   column of X contains the measured trait values for one species (tip in tree).
#'   Depending on the value of the global option "PCMBase.Internal.PC.Full" \code{NA}
#'   entries are interpreted either as missing measurements or non-existing traits
#'   (see \code{\link{PCMPresentCoordinates}}).
#' @param tree a phylo object with N tips.
#' @param model an S3 object specifying both, the model type (class, e.g. "OU") as
#'   well as the concrete model parameter values at which the likelihood is to be
#'   calculated (see also Details).
#' @param metaI a list returned from a call to \code{PCMValidate(tree, model)},
#'   containing meta-data such as N, M and k.
#' @param pruneI a named list containing cached preprocessing data for the tree used
#'   to perform post-order traversal (pruning). By default, this is created
#'   using \code{PCMPruningOrder(tree)}. This will use the default R-implementation of the
#'   likelihood calculation, which is based on the default R-implementation of the
#'   function \code{\link{PCMLmr}} (\code{\link{PCMLmr.default}}) and the S3 specification
#'   of the function
#'   \code{\link{PCMAbCdEf}} for the given model (a function called \code{PCMAbCdEf.MODEL},
#'   where model is the name of the mode and the class attribute of the \code{model}
#'   argument. For a different implementation of the function \code{\link{PCMLmr}},
#'   provide an S3 object for which an S3 specification has been written (see Details
#'   and example section).
#' @param pc a logical k x M matrix returned by \code{\link{PCMPresentCoordinates}}.
#' @param log logical indicating whether a log-liklehood should be calculated. Default
#'  is TRUE.
#' @param verbose logical indicating if some debug-messages should printed.
#'
#' @return a numerical value with named attributes as follows:
#' \item{X0}{A numerical vector of length k specifying the value at the root for which
#' the likelihood value was calculated. If the model contains a member called X0, this
#' vector is used; otherwise the value of X0 maximizing the likelihood for the given
#' model parameters is calculated by maximizing the quadratic polynomial
#' 'X0 * L_root * X0 + m_root * X0 + r_root'.}
#' \item{error}{A named list containing error information if a numerical or other
#' logical error occured during likelihood calculation (this is a list returned by
#'  \code{\link{PCMParseErrorMessage}}.}
#'  If an error occured during likelihood calculation, the default behavior is to
#'  return NA with a non-NULL error attribute. This behavior can be changed in
#'  using global options:
#'  \item{"PCMBase.Value.NA"}{Allows to specify a different NA value such as \code{-Inf} or \code{-1e20} which can be used in compbination with \code{log = TRUE} when
#'   using \code{optim} to maximize the log-likelihood;}
#'  \item{"PCMBase.Errors.As.Warnings"}{Setting this option to FALSE will cause any
#'  error to result in calling the \code{\link{stop}} R-base function. If not caught
#'  in a \code{\link{tryCatch}}, this will cause the inference procedure to abort at the occurence of a numerical error. By default, this option is set to TRUE, which
#'  means that \code{getOption("PCMBase.Value.NA", as.double(NA))} is returned with
#'  an error attribute and a warning is issued.}
#'
#' @details For efficiency, the arguments \code{metaI}, \code{pruneI} and \code{pc}
#'   can be provided explicitly, because these are not supposed to change during a
#'   model inference procedure such as likelihood maximization (see example).
#'
#' @seealso \code{\link{PCMValidate}} \code{\link{PCMPresentCoordinates}} \code{\link{PCMPruningOrder}} \code{\link{PCMAbCdEf}} \code{\link{PCMLmr}} \code{\link{PCMSim}} \code{\link{PCMCond}} \code{\link{PCMParseErrorMessage}}
#' @export
PCMLik <- function(
  X, tree, model,
  metaI = PCMInfo(X, tree, model, verbose = verbose),
  log = TRUE,
  verbose = FALSE) {

  # will change this value if there is no error
  value.NA <- getOption("PCMBase.Value.NA", as.double(NA))

  PCMLmr <- try(PCMLmr(X, tree, model, metaI, verbose = verbose, root.only = TRUE),
             silent = TRUE)

  if(class(PCMLmr) == "try-error") {
    errL <- PCMParseErrorMessage(PCMLmr)
    if(is.null(errL)) {
      err <- paste0("ERR:02041:PCMBase:MultivariatePCM.R:PCMLik:: There was a problem calculating the coefficients L,m,r. Error message from call to PCMLmr: ", PCMLmr, "; print(model):",
                    do.call(paste, c(as.list(capture.output(print(model))), list(sep="\n"))))
      errL <- PCMParseErrorMessage(err)
    } else {
      err <- PCMLmr
    }
    if(getOption("PCMBase.Errors.As.Warnings", TRUE)) {
      warning(err)
    } else {
      stop(err)
    }
    X0 <- model$X0
    attr(value.NA, 'X0') <- X0
    attr(value.NA, "error") <- errL

    return(value.NA)

  } else if(is.list(PCMLmr)) {

    L_root <- PCMLmr$L
    m_root <- PCMLmr$m
    r_root <- PCMLmr$r

    if(is.null(L_root) | is.null(m_root) | is.null(r_root)) {
      err <- paste0("ERR:02042:PCMBase:MultivariatePCM.R:PCMLik:: The list returned by PCMLmr did not contain elements 'L', 'm' and 'r'.")
      if(getOption("PCMBase.Errors.As.Warnings", TRUE)) {
        warning(err)
      } else {
        stop(err)
      }

      errL <- PCMParseErrorMessage(err)

      X0 <- model$X0
      attr(value.NA, 'X0') <- X0
      attr(value.NA, "error") <- errL

      return(value.NA)

    }

    if(is.null(model$X0)) {
      # set the root value to the one that maximizes the likelihood
      X0 <- try(solve(a=L_root + t(L_root), b = -m_root), silent = TRUE)
      if(class(X0) == "try-error") {
        err <- paste0(
          "ERR:02043:PCMBase:MultivariatePCM.R:PCMLik:: There was a problem calculating X0 from the coefficients L,m,r. ", "L=", toString(L), "; m=", toString(m), "; r=", r,
          ". Error message from call to solve(a=L_root + t(L_root), b = -m_root):", X0, "; print(model):",
          do.call(paste, c(as.list(capture.output(print(model))), list(sep="\n"))))

        errL <- PCMParseErrorMessage(err)
        if(getOption("PCMBase.Errors.As.Warnings", TRUE)) {
          warning(err)
        } else {
          stop(err)
        }
        X0 <- NULL
        attr(value.NA, "X0") <- X0
        attr(value.NA, "error") <- errL

        return(value.NA)

      }

    } else {

      X0 <- model$X0

    }

    loglik <- try(X0 %*% L_root %*% X0 + m_root %*% X0 + r_root, silent = TRUE)
    if(class(loglik) == "try-error") {
      err <- paste0(
        "ERR:02044:PCMBase:MultivariatePCM.R:PCMLik:: There was a problem calculating loglik from X0 and the coefficients L,m,r. ", "X0=", toString(X0), "L=", toString(L), "; m=", toString(m), "; r=", r,
        ". Error message from call to X0 %*% L_root %*% X0 + m_root %*% X0 + r_root:", loglik, "; print(model):",
        do.call(paste, c(as.list(capture.output(print(model))), list(sep="\n"))))

      errL <- PCMParseErrorMessage(err)
      if(getOption("PCMBase.Errors.As.Warnings", TRUE)) {
        warning(err)
      } else {
        stop(err)
      }

      attr(value.NA, "X0") <- X0
      attr(value.NA, "error") <- errL

      return(value.NA)

    }

    value <- try(as.vector(if(log) loglik else exp(loglik)), silent = TRUE)

    if(class(value) == "try-error") {
      err <- paste0(
        "ERR:02045:PCMBase:MultivariatePCM.R:PCMLik:: There was a problem calculating value from loglik=", toString(loglik), ". Error message from call to as.vector(if(log) loglik else exp(loglik)):", value, "; print(model):",
        do.call(paste, c(as.list(capture.output(print(model))), list(sep="\n"))))

      errL <- PCMParseErrorMessage(err)
      if(getOption("PCMBase.Errors.As.Warnings", TRUE)) {
        warning(err)
      } else {
        stop(err)
      }

      attr(value.NA, "X0") <- X0
      attr(value.NA, "error") <- errL

      return(value.NA)

    } else if(is.na(value)) {

      err <- paste0(
        "ERR:02046:PCMBase:MultivariatePCM.R:PCMLik:: There was a possible numerical problem, e.g. division of 0 by 0 when calculating the likelihood. value=", toString(value), "; calculated loglik=", toString(loglik), "; print(model):",
        do.call(paste, c(as.list(capture.output(print(model))), list(sep="\n"))), ". No error message was returned from the call to PCMLmr. Check for runtime warnings.")

      errL <- PCMParseErrorMessage(err)
      if(getOption("PCMBase.Errors.As.Warnings", TRUE)) {
        warning(err)
      } else {
        stop(err)
      }

      attr(value.NA, "X0") <- X0
      attr(value.NA, "error") <- errL

      return(value.NA)

    }

    # no errors detected, returning the calculated likelihood value:
    attr(value, "X0") <- X0
    return(value)
  }
}


#' Breadth-first tree traversal
#' @param tree a phylo object with possible singleton nodes (i.e. internal nodes with
#' one daughter node)
#' @return a vector of indices of edges in tree$edge in breadth-first order.
PCMPreorder <- function(tree) {
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
PCMPruningOrder <- function(tree) {
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


#' Determine which traits are present (active) on each node of the tree
#'
#' @description For every node (root, internal or tip) in tree, build a logical
#' vector of length k with TRUE values for every present coordinate. Non-present
#' coordinates arize from NA-values in the trait data. These can occur in two cases:
#' \describe{
#' \item{Missing measurements for some traits at some tips:}{the present coordinates
#' are FALSE for the corresponding tip and trait, but are full for all traits
#' at all internal and root nodes.}
#' \item{non-existent traits for some species:}{the FALSE present coordinates
#' propagate towards the parent
#' nodes - an internal or root node will have a present coordinate set to FALSE
#' for a given trait, if all of its descendants have this coordinate set to FALSE.}
#' }
#' These two cases have different effect on the likelihood calculation: missing
#' measurements are integrated out at the parent nodes; while non-existent traits
#' are treated as reduced dimensionality of the vector at the parent node.
#' The PCMBase package allows to specify the treatment of NA trait values using the
#' global option "PCMBase.Internal.PC.Full", which is set to TRUE by default. This
#' setting corresponds to assuming that all traits exist for all tips and NA values
#' are integrated out during likelihood calculation. Setting this option to FALSE
#' corresponds to assuming all NA values as non-existent.
#'
#' @param X numeric k x N matrix of observed values, with possible NA entries. The
#' columns in X are in the order of tree$tip.label
#' @param tree a phylo object
#' @param metaI a list returned by the PCMInfo function. Either leave this
#' as default or pass a previously computed metaI for the same X, tree and model.
#' @return a k x M logical matrix which can be passed as a pc argument to the PCMLik
#' function. The function fails in case when all traits are NAs for some of the tips.
#' In that case an error message is issued
#' "ERR:02001:PCMBase:MultivariatePCM.R:PCMPresentCoordinates:: Some tips have 0
#' present coordinates. Consider removing these tips.".
#' @seealso \code{\link{PCMLik}}
#' @export
PCMPresentCoordinates <- function(X, tree, metaI=PCMPruningOrder(tree)) {

  N <- metaI$N
  M <- metaI$M
  k <- metaI$k

  edge <- tree$edge
  endingAt <- metaI$endingAt
  nodesVector <- metaI$nodesVector
  nodesIndex <- metaI$nodesIndex
  nLevels <- metaI$nLevels
  unVector <- metaI$unVector
  unIndex <- metaI$unIndex
  unJ <- 1

  if(is.null(X)) {
    pc <- rep(TRUE, k*M)
  } else if(getOption("PCMBase.Internal.PC.Full", TRUE)) {
    pc <- rep(TRUE, k*M)
    dim(pc) <- c(k, M)
    pc[, 1:N] <- !is.na(X[, 1:N])
  } else {
    pc <- rep(FALSE, k*M)
    dim(pc) <- c(k, M)

    for(i in 1:nLevels) {
      nodes <- nodesVector[(nodesIndex[i]+1):nodesIndex[i+1]]
      es <- endingAt[nodes]

      if(nodes[1] <= N) {
        # all es pointing to tips
        pc[, nodes] <- !is.na(X[, nodes])
      } else {
        # edges pointing to internal nodes, for which all children nodes have been
        # visited
        # here we do nothing
      }

      #update parent pifs
      while(length(es)>0) {
        un <- unVector[(unIndex[unJ]+1):unIndex[unJ+1]]
        unJ <- unJ+1
        pc[, edge[es[un], 1]] <- pc[, edge[es[un], 1]] | pc[, edge[es[un], 2]]
        es <- es[-un]
      }
    }
    if(any(rowSums(pc) == 0)) {
      stop("ERR:02001:PCMBase:MultivariatePCM.R:PCMPresentCoordinates:: Some tips
           have 0 present coordinates. Consider removing these tips.")
    }
  }
  pc
}

#' Meta-information about a tree associated with a PCM
#'
#' @description
#' @inheritParams PCMLik
#' @return a named list with the following elements:
#' \item{M}{total number of nodes in the tree;}
#' \item{N}{number of tips;}
#' \item{k}{number of traits;}
#' \item{R}{number of regimes;}
#' \item{p}{number of free parameters describing the model;}
#' \item{r}{an integer vector corresponding to tree$edge with the regime for each
#' branch in tree;}
#' \item{xi}{an integer vector of 0's and 1's corresponding to the rows in tree$edge
#' indicating the presence of a jump at the corresponding branch;}
#' \item{pc}{a logical matrix of dimension k x M denoting the present coordinates
#' for each node;}
#'
#' @export
PCMInfo <- function(X, tree, model, verbose = FALSE) {
  UseMethod("PCMInfo", model)
}

#' Default PCMInfo implementation
#' @describeIn PCMInfo
#' @export
PCMInfo.PCM <- function(X, tree, model, verbose = FALSE) {
  res <- list(
    M = PCMNumNodes(tree),
    N = PCMNumTips(tree),
    k = PCMNumTraits(model),
    R = PCMNumRegimes(tree),
    p = PCMNumParams(model),
    r = PCMRegimes(tree, model),
    xi = PCMJumps(tree),
    preorder = PCMPreorder(tree)
  )
  res <- c(res, PCMPruningOrder(tree))
  res$pc <- PCMPresentCoordinates(X, tree, res)
  res
}

#' Multivariate normal distribution for a given tree and model
#' @description An S3 generic function that has to be implemented for every
#'  model class.
#' @inheritParams PCMLik
#' @param r an integer specifying a model regime
#' @return a list with the following members:
#' \item{omega}{}
#' \item{Phi}{}
#' \item{V}{}
#' @export
PCMCond <- function(tree, model, r=1, metaI=PCMInfo(NULL, tree, model, verbose), verbose = FALSE) {
  UseMethod("PCMCond", model)
}

#' Quadratic polynomial parameters A, b, C, d, E, f for each node
#' @description An S3 generic function that has to be implemented for every
#'  model class. This function is called by \code{\link{PCMLik}}.
#' @inheritParams PCMLik
#' @export
PCMAbCdEf <- function(tree, model, metaI=PCMInfo(NULL, tree, model, verbose), verbose = FALSE) {
  UseMethod("PCMAbCdEf", model)
}

#' Defaut R-implementation of PCMAbCdEf
#' @inheritParams PCMLik
#' @export
PCMAbCdEf.default <- function(tree, model, metaI=PCMInfo(NULL, tree, model, verbose), verbose = FALSE) {
  # number of regimes
  R <- metaI$R
  # number of tips
  N <- metaI$N
  # number of traits (variables)
  k <- metaI$k
  # number of nodes
  M <- metaI$M
  # present coordinates
  pc <- metaI$pc

  cond <- list()

  for(r in 1:R) {
    # create the conditional distribution for regime r
    cond[[r]] <- PCMCond(tree, model, r, metaI, verbose)
  }

  omega <- array(NA, dim=c(k, M))
  Phi <- array(NA, dim=c(k, k, M))
  V <- array(NA, dim=c(k, k, M))
  V_1 <- array(NA, dim=c(k, k, M))


  # returned general form parameters
  A <- array(NA, dim=c(k, k, M))
  b <- array(NA, dim=c(k, M))
  C <- array(NA, dim=c(k, k, M))
  d <- array(NA, dim=c(k, M))
  E <- array(NA, dim=c(k, k, M))
  f <- array(NA, dim=c(M))

  # vector of regime indices for each branch
  r <- metaI$r

  # identity k x k matrix
  I <- diag(k)

  # iterate over the edges
  for(edgeIndex in 1:(M-1)) {
    # parent node
    j <- tree$edge[edgeIndex, 1]
    # daughter node
    i <- tree$edge[edgeIndex, 2]

    # length of edge leading to i
    ti <- tree$edge.length[edgeIndex]

    # present coordinates in parent and daughte nodes
    kj <- pc[,j]
    ki <- pc[,i]

    omega[,i] <- cond[[r[edgeIndex]]]$omega(ti, edgeIndex)
    Phi[,,i] <- cond[[r[edgeIndex]]]$Phi(ti, edgeIndex)
    V[,,i] <- cond[[r[edgeIndex]]]$V(ti, edgeIndex)

    if(i<=N) {
      # add environmental variance at each tip node
      V[,,i] <- V[,,i] + model$Sigmae[,,r[edgeIndex]]
    }

    V_1[ki,ki,i] <- solve(V[ki,ki,i])

    `%op%` <- if(sum(ki) > 1) `%*%` else `*`

    A[ki,ki,i] <- (-0.5*V_1[ki,ki,i])
    E[kj,ki,i] <- t(Phi[ki,kj,i]) %op% V_1[ki,ki,i]
    b[ki,i] <- V_1[ki,ki,i] %*% omega[ki,i]
    C[kj,kj,i] <- -0.5 * matrix(E[kj,ki,i], sum(kj), sum(ki)) %*% matrix(Phi[ki,kj,i], sum(ki), sum(kj))
    d[kj,i] <- -E[kj,ki,i] %op% omega[ki,i]
    f[i] <- -0.5 * (t(omega[ki,i]) %*% V_1[ki,ki,i] %*% omega[ki,i] + sum(ki)*log(2*pi) + log(det(as.matrix(V[ki,ki,i]))))
  }

  list(A=A, b=b, C=C, d=d, E=E, f=f, omega=omega, Phi=Phi, V=V, V_1=V_1)
}

#' Quadratic polynomial parameters L, m, r
#'
#' @description
#'
#' @inheritParams PCMLik
#' @param root.only logical indicatin whether to return the calculated values of L,m,r
#'  only for the root or for all nodes in the tree.
#' @return A list with the members A,b,C,d,E,f,L,m,r for all nodes in the tree or
#'   only for the root if root.only=TRUE.
#' @export
PCMLmr <- function(X, tree, model,
                metaI = PCMInfo(X, tree, model, verbose = verbose),
                root.only = TRUE, verbose = FALSE) {
  UseMethod("PCMLmr", metaI)
}


#' Old defaut R-implementation of PCMLmr
#' @details This funciton is not a generic S3 implementation of PCMLmr - it needs
#' additional raguments and is designed to be called by PCMLik. It can still
#' be called by the end-user for debugging purpose.
#' @inheritParams PCMLik
#' @return A list with the members A,b,C,d,E,f,L,m,r for all nodes in the tree or
#'   only for the root if root.only=TRUE.
#' @export
PCMLmr.default <- function(
  X, tree, model,
  metaI = PCMInfo(X, tree, model, verbose = verbose),
  root.only = FALSE,
  verbose = FALSE) {

  unJ <- 1

  N <- metaI$N; M <- metaI$M; k <- metaI$k;

  edge <- tree$edge
  endingAt <- metaI$endingAt
  nodesVector <- metaI$nodesVector
  nodesIndex <- metaI$nodesIndex
  nLevels <- metaI$nLevels
  unVector <- metaI$unVector
  unIndex <- metaI$unIndex
  pc <- metaI$pc

  threshold_SV <- getOption("PCMBase.Threshold.SV", 1e-6)

  L <- array(0, dim=c(k, k, M))
  m <- array(0, dim=c(k, M))
  r <- array(0, dim=c(M))


  PCMAbCdEf <- PCMAbCdEf(tree = tree, model = model, metaI = metaI, verbose = verbose)


  # avoid redundant calculation
  log2pi <- log(2*pi)

  for(level in 1:nLevels) {
    nodes <- nodesVector[(nodesIndex[level]+1):nodesIndex[level+1]]
    es <- endingAt[nodes]

    if(nodes[1] <= N) {
      # all es pointing to tips
      #L[nodes,,] <- PCMAbCdEf$C[nodes,,]

      for(edgeIndex in es) {
        # parent and daughter nodes
        j <- edge[edgeIndex, 1]; i <- edge[edgeIndex, 2];
        # present coordinates
        kj <- pc[, j]; ki <- pc[, i];


        # check that V[ki,ki,] is non-singular
        svdV = svd(matrix(PCMAbCdEf$V[ki,ki,i], sum(ki)), 0, 0)$d
        if(min(svdV)/max(svdV) < threshold_SV) {
          err <- paste0(
            "ERR:02031:PCMBase:MultivariatePCM.R:PCMLmr.default:",i,":",
            " The matrix V for node ", i,
            " is nearly singular: min(svdV)/max(svdV)=", min(svdV)/max(svdV),
            ". Check the model parameters and the length of the branch",
            " leading to the node. For details on this error, read the User Guide.")
          stop(err)
        }

        # ensure symmetry of L[,,i]
        L[,,i] <- 0.5 * (PCMAbCdEf$C[,,i] + t(PCMAbCdEf$C[,,i]))

        r[i] <- with(PCMAbCdEf, t(X[ki,i]) %*% A[ki,ki,i] %*% X[ki,i] +
                       t(X[ki,i]) %*% b[ki,i] + f[i])

        m[kj,i] <- with(PCMAbCdEf, d[kj,i] + matrix(E[kj,ki,i], sum(kj), sum(ki)) %*% X[ki,i])

        #logDetV[i] <- with(PCMAbCdEf, det(-2*(A[i,ki,ki]+L[i,ki,ki])))

        #K <- K + sum(ki)
      }
    } else {
      # edges pointing to internal nodes, for which all children
      # nodes have been visited
      for(edgeIndex in es) {
        # parent and daughter nodes
        j <- edge[edgeIndex, 1]; i <- edge[edgeIndex, 2];
        # present coordinates
        kj <- pc[, j]; ki <- pc[, i];

        # check that V[i,ki,ki] is non-singular
        svdV = svd(matrix(PCMAbCdEf$V[ki,ki,i], sum(ki)), 0, 0)$d
        if(min(svdV)/max(svdV) < threshold_SV) {
          err <- paste0(
            "ERR:02031:PCMBase:MultivariatePCM.R:PCMLmr.default:",i,":",
            " The matrix V for node ", i,
            " is nearly singular: min(svdV)/max(svdV)=", min(svdV)/max(svdV),
            ", det(V)=", det(matrix(PCMAbCdEf$V[i,ki,ki], sum(ki))),
            ". Check the model parameters and the length of the branch",
            " leading to the node. For details on this error, read the User Guide.")
          stop(err)
        }



        # auxilary variables to avoid redundant evaluation
        AplusL <- as.matrix(PCMAbCdEf$A[ki,ki,i] + L[ki,ki,i])
        # ensure symmetry of AplusL, this should guarantee that AplusL_1 is symmetric
        # as well (unless solve-implementation is buggy.)
        AplusL <- 0.5 * (AplusL + t(AplusL))

        AplusL_1 <- solve(AplusL)

        EAplusL_1 <- matrix(PCMAbCdEf$E[kj,ki,i], sum(kj), sum(ki)) %*% AplusL_1
        logDetVNode <- log(det(-2*AplusL))

        # here it is important that we first evaluate r[i] and then m[i,kj]
        # since the expression for r[i] refers to to the value of m[i,ki]
        # before updating it.
        r[i] <- with(PCMAbCdEf, f[i]+r[i]+(sum(ki)/2)*log2pi-.5*logDetVNode -
                       .25*t(b[ki,i]+m[ki,i]) %*% AplusL_1 %*% (b[ki,i]+m[ki,i]))

        m[kj,i] <- with(PCMAbCdEf, d[kj,i] - .5*EAplusL_1 %*% (b[ki,i]+m[ki,i]))

        L[kj,kj,i] <- with(
          PCMAbCdEf,
          C[kj,kj,i] -.25*EAplusL_1 %*% t(matrix(E[kj,ki,i], sum(kj), sum(ki))))

        # ensure symmetry of L:
        L[kj,kj,i] <- 0.5 * (L[kj,kj,i] + t(L[kj,kj,i]))
      }
    }

    # add up to parents
    while(length(es)>0) {
      un <- unVector[(unIndex[unJ]+1):unIndex[unJ+1]]
      unJ <- unJ+1
      L[,,edge[es[un], 1]] <- L[,,edge[es[un], 1]] + L[,,edge[es[un], 2]]
      m[,edge[es[un], 1]] <- m[,edge[es[un], 1]] + m[,edge[es[un], 2]]
      r[edge[es[un], 1]] <- r[edge[es[un], 1]] + r[edge[es[un], 2]]
      es <- es[-un]
    }
  }

  if(root.only) {
    list(L = L[,,N+1],
         m = m[,N+1],
         r = r[N+1])
  } else {
    c(PCMAbCdEf[c("A", "b", "C", "d", "E", "f", "V", "V_1")],
      list(L = L, m = m, r = r))
  }
}



#' Extract error information from a formatted error message.
#' @param x character string representing the error message.
#' @description The function searches x for a pattern matching the format
#' 'ERR:5-alphanumeric-character-code:project-name:source-file:error-specifics:'.
#' Specifically it
#' searches for a regular expression pattern "ERR:[0-9a-zA-Z]+:[^:]+:[^:]+:[^:]+:[^:]*:".
#' @return a named list with the parsed error information or NULL, if no match
#' was found. The elements of this list are named as follows:
#' \item{type}{The type of the error message. Usually this is ERROR, but could be
#' WARNING or anything else.}
#' \item{icode}{An an alphanumeric code of the error.}
#' \item{project}{The name of the project locating the code that raised the error.}
#' \item{file}{The name of the source-file containing the code that raised the error.}
#' \item{fun}{The name of the function raising the error}
#' \item{info}{A character string containing additional error-specific information}
#' \item{msg}{A verbal description of the error.}
#' @export
PCMParseErrorMessage <- function(x) {
  res <- try({
    if(is.character(x)) {
      code <- regmatches(x, regexpr("ERR:[0-9a-zA-Z]+:[^:]+:[^:]+:[^:]+:[^:]*:", x))
      if(length(code) > 0) {
        code <- code[1]
        codeL <- strsplit(code, split=":")[[1]]
        list(
          type = codeL[1],
          icode = codeL[2],
          project = codeL[3],
          file = codeL[4],
          fun = codeL[5],
          info = codeL[6],
          code = code,
          msg = x
        )
      } else {
        NULL
      }
    } else {
      NULL
    }
  }, silent = TRUE)

  if(class(res)=="try-error") {
    NULL
  } else {
    res
  }
}


#' Sums of pairs of elements in a vector
#' @param lambda a numeric vector
#' @return a squared symmetric matrix with elem_ij=lambda_i+lambda_j.
#'
PCMPairSums <- function(lambda) {
  sapply(lambda, function(li) sapply(lambda, function(lj) li+lj))
}

#' Eigen-decomposition of a matrix H
#' @param H a numeric matrix
#' @return a list with elements as follows:
#' \item{lambda}{a vector of the eigenvalues of H}
#' \item{P}{a squared matrix with column vectors, the eigenvectors of H corresponding to the
#' eigenvalues in \code{lambda}}
#' \item{P_1}{the inverse matrix of \code{P}}.
#' @details The function fails with an error message if H is defective, that is, if its matrix of
#' eigenvectors is computationally singular. The test for singularity is based on the \code{\link{rcond}} function.
#'
PCMPLambdaP_1 <- function(H) {
  # here the argument H is a matrix specifying the alphas in a OU process
  r <- eigen(H)
  if(isTRUE(all.equal(rcond(r$vectors),0))) {
    stop(paste0("ERR:02041:PCMBase:MultivariatePCM.R:PCMPLambdaP_1:: The provided H matrix is defective. Its matrix of eigenvectors is computationally singular. H=", toString(H)))
  }
  list(lambda=r$values, P=r$vectors, P_1=solve(r$vectors))
}


#' Create a function of time that calculates (1-exp(-lambda_{ij}*time))/lambda_{ij}
#' for every element lambda_{ij} of the input matrix Lambda_ij.
#'
#' @param Lambda_ij a squared numerical matrix of dimension k x k
#' @param threshold.Lambda_ij a 0-threshold for abs(Lambda_i + Lambda_j), where Lambda_i
#'   and Lambda_j are eigenvalues of the parameter matrix H. This
#'   threshold-value is used as a condition to
#'   take the limit time of the expression `(1-exp(-Lambda_ij*time))/Lambda_ij` as
#'   `(Lambda_i+Lambda_j) --> 0`. You can control this value by the global option
#'   "PCMBase.Threshold.Lambda_ij". The default value (1e-8) is suitable for branch
#'    lengths bigger than 1e-6. For smaller branch lengths, you are may want to
#'    increase the threshold value using, e.g.
#'   `options(PCMBase.Threshold.Lambda_ij=1e-6)`.
#' @details the function (1-exp(-lambda_{ij}*time))/lambda_{ij} corresponds to the product
#' of the CDF of an exponential distribution with rate Lambda_{ij} multiplied by its mean value (mean waiting time).
#' @return a function of time returning a matrix with entries formed from the
#'  above function or the limit, time, if |Lambda_{ij}|<=trehshold0.
PCMPExpxMeanExp <- function(
  Lambda_ij,
  threshold.Lambda_ij = getOption("PCMBase.Threshold.Lambda_ij", 1e-8) ) {

  idx0 <- which(abs(Lambda_ij)<=threshold.Lambda_ij)
  function(time) {
    res <- (1-exp(-Lambda_ij*time))/Lambda_ij
    if(length(idx0)>0) {
      res[idx0] <- time
    }
    res
  }
}


#' Variance-covariance matrix of an OU process with optional jump at the start
#' @param H a numerical k x k matrix - selection strength parameter.
#' @param Sigma a numerical k x k matrix - neutral dift unit-time variance-covariance matrix.
#' @param Sigmaj is the variance matrix of the normal jump distribution (default is NULL).
#' @param xi a vector of 0's and 1's corresponding to each branch in the tree. A value of 1
#' indicates that a jump takes place at the beginning of the branch. This arugment is only
#' used if Sigmaj is not NULL. Default is NULL.
#'
#' @param threshold.Lambda_ij a 0-threshold for abs(Lambda_i + Lambda_j), where Lambda_i
#' and Lambda_j are eigenvalues of the parameter matrix H. This threshold-values is used as
#' a condition to take the limit time of the expression `(1-exp(-Lambda_ij*time))/Lambda_ij`
#' as `(Lambda_i+Lambda_j) --> 0`. You can control this value by the global option
#' "PCMBase.Threshold.Lambda_ij". The default value (1e-8) is suitable for branch lengths
#' bigger than 1e-6. For smaller branch lengths, you may want to increase the threshold
#' value using, e.g.  `options(PCMBase.Threshold.Lambda_ij=1e-6)`.
#' @return a function of one numerical argument (time) and an integer indicating the branch-index
#' that is used to check the corresponding element in xi.
#' @export
PCMCondVOU <- function(
  H, Sigma , Sigmaj=NULL, xi=NULL,
  e_Ht = NULL,
  threshold.Lambda_ij = getOption("PCMBase.Threshold.Lambda_ij", 1e-8)) {

  force(H)
  force(Sigma)
  force(Sigmaj)
  force(xi)
  force(e_Ht)

  if(is.null(dim(H)) | is.null(dim(Sigma))) {
    stop('ERR:02051:PCMBase:MultivariatePCM.R:PCMCondVOU:: H and Sigma must be k x k matrices.')
  }

  k <- dim(Sigma)[1]

  if(!is.matrix(Sigma) | !is.matrix(H) | !isTRUE(all.equal(c(dim(H), dim(Sigma)), rep(k, 4)))) {
    stop(paste0("ERR:02052:PCMBase:MultivariatePCM.R:PCMCondVOU:: H and Sigma must be  ", k, " x ", k, " matrices."))
  }

  if(!is.null(Sigmaj) & !is.matrix(Sigmaj) & !isTRUE(all.equal(dim(Sigmaj), rep(k, 2)))) {
    stop(paste0("ERR:02053:PCMBase:MultivariatePCM.R:PCMCondVOU:: Sigmaj must be NULL or a ", k, " x ", k, " matrix."))
  }

  if(!is.null(e_Ht) & !is.matrix(e_Ht) & !isTRUE(all.equal(dim(e_Ht), rep(k, 2)))) {
    stop(paste0("ERR:02054:PCMBase:MultivariatePCM.R:PCMCondVOU:: e_Ht must be NULL or a ", k, " x ", k, " matrix."))
  }

  PLP_1 <- PCMPLambdaP_1(H)

  Lambda_ij <- PCMPairSums(PLP_1$lambda)
  fLambda_ij <- PCMPExpxMeanExp(Lambda_ij, threshold.Lambda_ij)
  P_1SigmaP_t <- PLP_1$P_1 %*% Sigma %*% t(PLP_1$P_1)


  function(t, edgeIndex, e_Ht = NULL) {
    res <- PLP_1$P %*% (fLambda_ij(t) * P_1SigmaP_t) %*% t(PLP_1$P)
    if(!is.null(Sigmaj)) {
      if(is.null(e_Ht)) {
        e_Ht <- expm(-t*H)
      }
      res <- res + xi[edgeIndex]*(e_Ht %*% Sigmaj %*% t(e_Ht))
    }
    Re(res)
  }
}

PCMCondRandom <- function(PCMCondObject, n=1, x0, t, edgeIndex) {
  with(PCMCondObject, {
    rmvnorm(n=n, mean = omega(t, edgeIndex) + Phi(t, edgeIndex)%*%x0, sigma=V(t, edgeIndex))
  })
}

PCMCondDensity <- function(PCMCondObject, x, x0, t, edgeIndex, log=FALSE) {
  with(PCMCondObject, {
    dmvnorm(x, mean=omega(t, edgeIndex) + Phi(t, edgeIndex)%*%x0, sigma=V(t, edgeIndex), log=log)
  })
}
