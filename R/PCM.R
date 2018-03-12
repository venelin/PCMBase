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
#' @title Create a phylogenetic comparative model object
#'
#' @description This is the entry-point function for creating model objects within
#' the PCMBase framework. A PCM
#' @param model This argument can take one of the following forms:
#' \itemize{
#' \item a character vector of the S3-classes of the model object to be
#' created (one model object can have one or more S3-classes, with the class
#' PCM at the origin of the hierarchy);
#' \item an S3 object which's class inherits from the PCM S3 class.
#' }
#' The Details section explains how these two types of input are processed.
#' @param k integer denoting the number of traits (defaults to 1).
#' @param regimes a character or integer vector denoting the regimes.
#' @param params NULL (default) or a list of parameter values (scalars, vectors,
#' matrices, or arrays) or sub-models (S3 objects inheriting from the PCM class).
#' See details.
#' @param vecParams NULL (default) or a numeric vector the vector
#' representation of the variable parameters in the model. See details.
#' @param offset integer offset in vecParams; see Details.
#' @param ... additional parameters intended for use by sub-classes of the PCM
#' class.
#'
#' @details
#' @return an object of S3 class as defined by the argument model.
#' @export
PCM <- function(model, k = 1, regimes = 1,
                params = NULL, vecParams = NULL, offset = 0, ...) {
  UseMethod("PCM", model)
}

#' @export
PCM.default <- function(model, k = 1, regimes = 1,
                        params = NULL, vecParams = NULL, offset = 0, ...) {
  stop(paste0("ERR:02091:PCMBase:PCM.R:PCM.default:: You should provide a PCM object, but provided a ", class(model)[1]))
}

#' Common constructor for all PCM classes
#' @inheritParams PCM
#' @export
PCM.PCM <- function(model, k = 1, regimes = 1,
                    params = NULL, vecParams = NULL, offset = 0, ...) {

  if(is.null(attr(model, "specParams", exact = TRUE))) {
    attr(model, "specParams") <- PCMSpecifyParams(model, ...)
  }

  if(is.null(attr(model, "p", exact = TRUE))) {
    attr(model, "p") <- PCMNumParams(model)
  }

  specParams <- attr(model, "specParams")

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

#' @inherit PCM
#' @export
PCM.character <- function(model, k = 1, regimes = 1,
                          params = NULL, vecParams = NULL, offset = 0, ...) {
  modelObj <- list()
  class(modelObj) <- model
  if(length(model) == 1) {
    class(modelObj) <- c(class(modelObj), PCMParentClasses(modelObj))
  }
  attr(modelObj, "k") <- k
  attr(modelObj, "regimes") <- regimes

  PCM(modelObj, k, regimes, params, vecParams, offset, ...)
}

#' @export
is.PCM <- function(x) inherits(x, "PCM")

#' @export
print.PCM <- function(x, ...) cat (format(x, ...), "\n")

#' @export
format.PCM <- function(x, ...) {
  if( !is.PCM(x) ) {
    stop("ERR:02061:PCMBase:PCM.R:format.PCM:: x must inherit from S3 class PCM.")
  }
  specParams <- attr(x, "specParams", exact = TRUE)

  res <- paste0(
    PCMDescribe(x, ...),
    "\nS3 class: ", toString(class(x)), "; ",
    "k=", attr(x, "k", exact = TRUE), "; p=", PCMNumParams(x, ...), "; ",
    "regimes: ", toString(attr(x, "regimes")), ". Parameters/sub-models:\n")

  for(name in names(x)) {
    if(!name %in% names(specParams)) {
      stop(paste0("ERR:02062:PCMBase:PCM.R:format.PCM:: model member name '",
                  name, "' not found among the names in model specParams attribute."))
    }

    type <- specParams[[name]]$type
    description <- specParams[[name]]$description

    strList <- as.list(capture.output(print(x[[name]])))
    strList$sep = "\n"
    res <- paste0(res, name, " (", toString(type),
                  if(!is.null(description)) paste0("; ", description) else "", "):\n",
                  do.call(paste, strList))
    res <- paste0(res, "\n")
  }
  res
}

#' Parent S3 classes for a model class
#' @param model an S3 object.
#' @details This S3 generic function is intended to be specified for user models.
#' This function is called by PCM.character to determine the parent classes for
#' a given model class.
#' @return a vector of character string denoting the names of the parent classes
#' @seealso \code{\link{PCM.character}}
#' @export
PCMParentClasses <- function(model) {
  UseMethod("PCMParentClasses", model)
}

#' Human friendly description of a PCM
#' @param model a PCM model object
#' @details This S3 generic function is intended to be specified for user models
#' @return a character string
#' @export
PCMDescribe <- function(model, ...) {
  UseMethod("PCMDescribe", model)
}

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

#' @export
PCMNumTraits.PCM <- function(model) {
  attr(model, "k")
}

#' Number of regimes in a model
#' @param model a PCM object
#' @return an integer
#' @export
PCMNumRegimes <- function(model) {
  UseMethod("PCMNumRegimes", model)
}

#' @export
PCMNumRegimes.PCM <- function(model) {
  length(attr(model, "regimes"))
}

#' Regimes in a model
#' @param model a PCM object
#' @return a character or an integer vector giving the regime names of the models
#' @export
PCMRegimes <- function(model) {
  UseMethod("PCMRegimes", model)
}

#' @export
PCMRegimes.PCM <- function(model) {
  attr(model, "regimes")
}

#' Number of free parameters describing fully a PCM
#' @param model a PCM object
#' @param ... other arguments (possible future use)
#' @return an integer
#' @export
PCMNumParams <- function(model, ...) {
  UseMethod("PCMNumParams", model)
}

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
      p <- p + PCMNumParams(specParams[[name]]$default, ...)
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
PCMSetDefaultRegime <- function(tree, regime) {
  if(is.PCM(regime)) {
    regime <- PCMRegimes(regime)
  }
  eval(substitute(tree$edge.regime <- rep(regime[1], length(tree$edge.length))), parent.frame())
}

#' Unique regimes on a tree
#' @description This is a shortcut for sort(unique(tree$edge.regime))
#' @param tree a phylo object with an additional member edge.regime which should
#' be a character or an integer vector of length equal to the number of branches.
#'
#' @return a character or an integer vector depending on tree$edge.regime.
#' @export
PCMUniqueRegimesTree <- function(tree) {
  if(is.null(tree$edge.regime)) {
    stop("ERR:02010:PCMBase:PCM.R:PCMRegimesUniqueTree:: tree$edge.regime is NULL,
         but should be a character or an integer vector denoting regime names.")
  }
  sort(unique(tree$edge.regime))
}

#' Number of unique regimes on a tree
#' @param tree a phylo object
#' @return the number of different regimes encountered on the tree branches
#' @export
PCMNumUniqueRegimesTree <- function(tree) {
  length(PCMUniqueRegimesTree(tree))
}

#' Jumps in modeled traits associated with branches in a tree
#' @inheritParams PCMNumTips
#' @return an integer vector of 0's and 1's with entries correspondin to the
#' denoting if a jump took place at the beginning of a branch.
PCMJumps <- function(tree) {
  if(!is.null(tree$edge.jump)) {
    if(length(tree$edge.jump) != nrow(tree$edge) |
       !isTRUE(all(tree$edge.jump %in% c(0L, 1L)))) {
      stop("ERR:02081:PCMBase:PCM.R:PCMJumps:: tree$edge.jump should
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
PCMMatchRegimesModelTree <- function(tree, model) {
  regimes <- match(tree$edge.regime, PCMRegimes(model))
  if(any(is.na(regimes))) {
    stop(paste0("ERR:02071:PCMBase:PCM.R:PCMMatchRegimesModelTree:: ",
                " Some of the regimes in tree$edge.regime not found in",
                "attr(model, 'regimes').\n",
                "unique regimes on the tree:", toString(PCMUniqueRegimesTree(tree)), "\n",
                "attr(model, 'regimes')", toString(PCMRegimes(model))))
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
#' \describe{
#' \item{"rep"}{repeat \code{vecParams[offset+1]} k times}
#' \item{"full"}{copy \code{vecParams[offset + (1:k)]}}
#' \item{"custom"}{use a custom template vector which elements specified by mask are
#' to be replaced with \code{vecParams[indices(offset, k)]}}
#' }.
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
      stop("ERR:020a1:PCMBase:PCM.R:PCMLoadVectorParameter:: indices should be a
           function(offset, k) returning an integer vector.")
    }
  } else if(type[1] == "fixed") {
    num <- 0
  } else {
    stop(paste0("ERR:020a2:PCMBase:PCM.R:PCMLoadVectorParameter:: type ", type[1], " not recognized."))
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
#' \describe{
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
#' }
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
      stop("ERR:020b1:PCMBase:PCM.R:PCMLoadMatrixParameter:: indices should be a
           function(offset, k) returning an integer vector.")
    }
  } else if(type[1] == "fixed") {
    num <- 0
  } else {
    stop(paste0("ERR:020b2:PCMBase:PCM.R:PCMLoadMatrixParameter:: type ", type[1], " not recognized."))
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

#' @export
PCMSetParams.PCM <- function(model, params, inplace = TRUE, ...) {
  specParams <- attr(model, "specParams")
  for(name in names(params)) {
    if(! (name%in%names(specParams)) ) {
      stop(paste0("ERR:020c1:PCMBase:PCM.R:PCMSetParams:: ", name,
                  " is not a settable parameter of the model, check attr(model, 'specParams')."))
    }
    type <- specParams[[name]]$type

    if(type[1] != "model") {
      if(! identical(length(model[[name]]), length(params[[name]])) ) {
        stop(paste0("ERR:020c2:PCMBase:PCM.R:PCMSetParams:: params[[", name,
                    "]] is not the same length as model[[", name, "]]; ",
                    "length(params[[", name, "]])=", length(params[[name]]),
                    ", length(model[[", name, "]])=", length(model[[name]]), "."))
      }
      if(! identical(dim(model[[name]]), dim(params[[name]])) ) {
        stop(paste0("ERR:020c3:PCMBase:PCM.R:PCMSetParams:: params[[", name,
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

#' Conditional distribution of a daughter node given its parent node
#' @description An S3 generic function that has to be implemented for every
#'  model class.
#' @inheritParams PCMLik
#' @param r an integer specifying a model regime
#' @return an object of type specific to the type of model
#' @export
PCMCond <- function(tree, model, r=1, metaI=PCMInfo(NULL, tree, model, verbose), verbose = FALSE) {
  UseMethod("PCMCond", model)
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
#' \describe{
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
#' "ERR:02002:PCMBase:PCM.R:PCMSim:: X0 must be of length ...".
#'}
#' @importFrom mvtnorm rmvnorm
#' @seealso \code{\link{PCMLik}} \code{\link{PCMInfo}} \code{\link{PCMCond}}
#' @export
PCMSim <- function(
  tree, model, X0,
  metaI = PCMInfo(X = NULL, tree = tree, model = model, verbose = verbose),
  verbose = FALSE) {
  UseMethod("PCMSim", model)
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
#' @param metaI a list returned from a call to \code{PCMInfo(X, tree, model)},
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
#' @param log logical indicating whether a log-liklehood should be calculated. Default
#'  is TRUE.
#' @param verbose logical indicating if some debug-messages should printed.
#'
#' @return a numerical value with named attributes as follows:
#' \describe{
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
#'}
#' @details For efficiency, the argument \code{metaI}
#'   can be provided explicitly, because this is not supposed to change during a
#'   model inference procedure such as likelihood maximization.
#'
#' @seealso \code{\link{PCMInfo}} \code{\link{PCMAbCdEf}} \code{\link{PCMLmr}} \code{\link{PCMSim}} \code{\link{PCMCond}} \code{\link{PCMParseErrorMessage}}
#' @export
PCMLik <- function(
  X, tree, model,
  metaI = PCMInfo(X, tree, model, verbose = verbose),
  log = TRUE,
  verbose = FALSE) {
  UseMethod("PCMLik", model)
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
#' "ERR:02001:PCMBase:PCM.R:PCMPresentCoordinates:: Some tips have 0
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
      stop("ERR:02001:PCMBase:PCM.R:PCMPresentCoordinates:: Some tips
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
#' \item{RTree}{number of regimes on the tree (distinct elements of tree$edge.regime);}
#' \item{RModel}{number of regimes in the model (distinct elements of attr(model, regimes));}
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

#' @export
PCMInfo.PCM <- function(X, tree, model, verbose = FALSE) {
  res <- list(
    M = PCMNumNodes(tree),
    N = PCMNumTips(tree),
    k = PCMNumTraits(model),
    RTree = PCMNumUniqueRegimesTree(tree),
    p = PCMNumParams(model),
    r = PCMMatchRegimesModelTree(tree, model),
    xi = PCMJumps(tree),
    edge = tree$edge,
    edge.length = tree$edge.length,
    preorder = PCMPreorder(tree)
  )
  res <- c(res, PCMPruningOrder(tree))
  res$pc <- PCMPresentCoordinates(X, tree, res)
  res
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
