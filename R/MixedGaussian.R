# Copyright 2016 Venelin Mitov
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

#' @export
PCMModelTypes.MixedGaussian <- function(obj) {
  attr(obj, "modelTypes", exact = TRUE)
}



#' @export
PCMParentClasses.MixedGaussian <- function(model) {
  c("GaussianPCM", "PCM")
}

#' @export
PCMDescribe.MixedGaussian <- function(model, ...) {
  "Mixed Gaussian model"
}

#' @export
PCMCond.MixedGaussian <- function(
  tree, model, r=1,
  metaI = PCMInfo(NULL, tree, model, verbose = verbose), verbose=FALSE) {
  if(! is.MixedGaussian(model) ) {
    stop("MixedGaussian.R:PCMCond.MixedGaussian:: model should inherit from S3 class 'MixedGaussian'.")
  }
  if(!is.character(r)) {
    # integer regime number must be mapped to a character regime because
    # the model object may contain global parameter vectors and matrices apart
    # from models associated with regimes.
    r <- as.character(PCMRegimes(model)[r])
  }

  if(!is.null(model$Sigmae_x)) {
    Sigmae <- as.matrix(model$Sigmae_x) %*% t(as.matrix(model$Sigmae_x))
  } else {
    Sigmae <- NULL
  }

  if(is.null(Sigmae)) {
    # Sigmae_x is omitted in the parent model, so just return the result from
    # the sub-model
    PCMCond(tree, model[[r]], 1, metaI, verbose = verbose)
  } else {
    OmegaPhiV2 <- OmegaPhiV <- PCMCond(tree, model[[r]], 1, metaI, verbose = verbose)
    if(is.null(model[[r]]$Sigmae_x)) {
      # Sigmae_x is _Omitted in model[[r]]
      OmegaPhiV2$V <- function(t, edgeIndex, metaI, e_Ht = NULL) {
        res <- OmegaPhiV$V(t, edgeIndex, metaI, e_Ht)
        if(metaI$edge[edgeIndex,2] <= metaI$N) {
          res <- res + Sigmae
        }
        res
      }
    } else {
      # Sigmae_x is not _Omitted in model[[r]], so it overwrites the global one.
      OmegaPhiV2$V <- function(t, edgeIndex, metaI, e_Ht = NULL) {
        OmegaPhiV$V(t, edgeIndex, metaI, e_Ht)
      }
    }

    OmegaPhiV2
  }
}

#' @export
PCMParamCount.MixedGaussian <- function(
  o, countRegimeChanges = FALSE, countModelTypes = FALSE,  offset = 0, k = 1,
  R = 1, parentModel = NULL) {

  p <- NextMethod()
  if(countModelTypes) {
    if(length(PCMModelTypes(o)) > 1) {
      # + 1 has already been counted by the default PCM implementatioin (call
      # to NextMethod() above)
      p <- p + PCMNumRegimes(o) - 1
    }
  }
  p
}

#' @export
PCMExtractRegimes.MixedGaussian <- function(obj, regimes = seq_len(PCMNumRegimes(obj))) {
  regimes.obj <- PCMRegimes(obj)
  regimes.obj2 <- regimes.obj[regimes]

  names.obj2 <- setdiff(names(obj), setdiff(regimes.obj, regimes.obj2))
  obj2 <- obj[names.obj2]

  for(na in names(attributes(obj))) {
    if(na != "names") {
      attr(obj2, na) <- attr(obj, na)
    }
  }

  attr(obj2, "regimes") <- regimes.obj2
  attr(obj2, "mapping") <- attr(obj2, "mapping")[regimes.obj2]
  attr(obj2, "p") <- PCMParamCount(obj2)

  names.spec.obj2 <- setdiff(names(attr(obj, "spec")), setdiff(regimes.obj, regimes.obj2))
  attr(obj2, "spec") <- attr(obj2, "spec")[names.spec.obj2]

  class(obj2)[1L] <- class(attr(obj2, "spec"))[1L] <-
    paste0("MixedGaussian", "_", do.call(paste0, as.list(regimes.obj2)))

  obj2
}

#' Create a multi-regime Gaussian model (MixedGaussian)
#' @param k integer specifying the number of traits.
#' @param modelTypes,mapping These two arguments define the mapping between the
#' regimes in the model and actual types of models. For convenience, different
#' combinations are possible as explained below:
#' \itemize{
#' \item \code{modelTypes} is a (possibly named) character string vector. Each
#' such string denotes a mixed Gaussian regime model class, e.g. the result of calling
#' \code{MGPMDefaultModelTypes()}. In that case \code{mapping} can be either
#' an integer vector with values corresponding to indices in \code{modelTypes}
#' or a character string vector. If \code{mapping} is a character string vector,
#' first it is matched against \code{names(modelTypes)} and if the match fails
#' either because of \code{names(modelTypes)} being \code{NULL} or because some
#' of the entries in \code{mapping} are not present in \code{names(modelTypes)},
#' then an attempt is made to match \code{mapping} against \code{modelTypes},
#' i.e. it is assumed that \code{mapping} contains actual class names.
#' \item \code{modelTypes} is a (possibly named) list of PCM models of
#' \code{k} traits. In this case \code{mapping} can again be an integer vector
#' denoting indices in \code{modelTypes} or a character string vector denoting
#' names in \code{modelTypes}.
#' }
#' As a final note, \code{mapping} can also be named. In this case the names are
#' assumed to be the names of the different regimes in the model. If
#' \code{mapping} is not named, the regimes are named automatically as
#' \code{as.character(seq_len(mapping))}.  For example, if
#' \code{modelTypes = c("BM", "OU")} and
#' \code{mapping = c(a = 1, b = 1, c = 2, d = 1)} defines an MixedGaussian with
#' four different regimes named 'a', 'b', 'c', 'd', and  model-types
#' BM, BM, OU and BM, corresponding to each regime.
#' @param className a character string defining a valid S3 class name for the
#' resulting MixedGaussian object. If not specified, a className is generated
#' using the expression
#' \code{ paste0("MixedGaussian_", do.call(paste0, as.list(mapping)))}.
#'
#' @param X0 specification for the global vector X0 to be used by all
#' models in the MixedGaussian.
#'
#' @param ... specifications for other _Global parameters coming after X0.
#'
#' @param Sigmae_x sepcification of a _Global Sigmae_x parameter. This is used
#' by Submodels only if they have Sigmae_x _Omitted.
#'
#' @return an object of S3 class className inheriting from MixedGaussian,
#' GaussianPCM and PCM.
#'
#' @details If X0 is not NULL it has no sense to use model-types including X0 as
#' a parameter (e.g. use BM1 or BM3 instead of BM or BM2). Similarly if Sigmae_x
#' is not NULL there is no meaning in using model-types including Sigmae_x as a
#' parameter, (e.g. use OU2 or OU3 instead of OU or OU1).
#' @seealso \code{\link{PCMTreeGetPartNames}}
#' @seealso \code{\link{PCMModels}()}
#'
#' @export
MixedGaussian <- function(
  k,
  modelTypes,
  mapping,
  className = paste0("MixedGaussian_", do.call(paste0, as.list(mapping))),
  X0 = structure(0.0, class = c("VectorParameter", "_Global"),
                 description = "trait values at the root"),
  ...,
  Sigmae_x = structure(
    0.0,
    class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
    description = "Upper triangular factor of the non-phylogenetic variance-covariance")) {

  regimes <- if(is.null(names(mapping))) {
    as.character(seq_along(mapping))
  } else {
    names(mapping)
  }
  if(is.character(mapping)) {
    if(is.null(names(modelTypes)) || any(is.na(mapping2 <- match(mapping, names(modelTypes)))) ) {
      mapping2 <- match(mapping, modelTypes)
    }
    if(any(is.na(mapping2))) {
      stop(
        paste0(
          "MixedGaussian:: some of the models in mapping not found in modelTypes: ",
          "modelTypes = ", toString(modelTypes),
          ", mapping =", toString(mapping)))
    } else {
      mapping <- mapping2
    }
  }

  mappingModelRegime <- modelTypes[mapping]
  names(mappingModelRegime) <-regimes

  modelTypesAreStrings <- FALSE
  spec <- list(X0 = X0, ...)
  for(m in seq_along(mapping)) {
    subModel <- mappingModelRegime[[m]]
    if(is.character(subModel)) {
      modelTypesAreStrings <- TRUE
      subModel <- structure(0.0, class = subModel)
      class(subModel) <- c(class(subModel), PCMParentClasses(subModel))
    }
    spec[[regimes[m]]] <- PCMSpecify(subModel)
    attr(spec[[regimes[m]]], "k") <- as.integer(k)
    attr(spec[[regimes[m]]], "regimes") <- 1L
  }
  spec[["Sigmae_x"]] <- Sigmae_x

  attr(spec, "k") <- as.integer(k)
  class(spec) <- c(className, "MixedGaussian", "GaussianPCM", "PCM")

  if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')

  if(is.list(modelTypes)) {
    res <- PCM(
      model = class(spec), modelTypes = sapply(modelTypes, function(m) class(m)[1L]),
      k = k, regimes = regimes, spec = spec, params = mappingModelRegime)
  } else if(modelTypesAreStrings) {
    res <- PCM(model = class(spec), modelTypes, k = k, regimes = regimes, spec = spec)
  } else {
    stop("MixedGaussian:: Model types should either be a character vector or a list of PCM objects.")
  }

  attr(res, "mapping") <- mapping
  attr(res, "spec") <- spec
  res
}

#' Check if an object is a `MixedGaussian` PCM
#' @param x any object
#' @return TRUE if x inherits from the S3 class `MixedGaussian`, FALSE otherwise.
#'
#' @export
is.MixedGaussian <- function(x) {
  inherits(x, "MixedGaussian")
}

#' Convert a \code{GaussianPCM} model object to a \code{MixedGaussian} model object
#' @param o an R object: either a \code{GaussianPCM} or a \code{MixedGaussian}.
#' @param modelTypes NULL (the default) or a (possibly named) character string
#' vector. Each such string denotes a mixed Gaussian regime model class, e.g.
#' the result of calling \code{MGPMDefaultModelTypes()}. If specified, an
#' attempt is made to match the deduced Gaussian regime model type from \code{o}
#' with the elements of \code{modelTypes} and an error is raised if the match
#' fails. If the match succeeds the converted MixedGaussian object will have the
#' specified \code{modelTypes} parameter as an attribute \code{"modelTypes"}.
#' @return a \code{MixedGaussian} object.
#' @examples
#' mg <- as.MixedGaussian(PCMBaseTestObjects$model.ab.123.bSigmae_x)
#' stopifnot(
#'   PCMLik(
#'     X = PCMBaseTestObjects$traits.ab.123,
#'     PCMBaseTestObjects$tree.ab,
#'     PCMBaseTestObjects$model.ab.123.bSigmae_x) ==
#'   PCMLik(
#'     X = PCMBaseTestObjects$traits.ab.123,
#'     PCMBaseTestObjects$tree.ab,
#'     mg))
#'
#' @export
as.MixedGaussian <- function(o, modelTypes = NULL) {
  if(is.MixedGaussian(o)) {
    if(!is.null(modelTypes)) {
      modelTypeso <- PCMModelTypes(o)
      m <- match(modelTypeso, modelTypes)
      if(isTRUE(any(is.na(m)))) {
        stop(paste0("as.MixedGaussian:: some of the modelTypes in o could not",
                    "be matched against the supplied modelTypes."))
      }
      if(!identical(modelTypeso, modelTypes)) {
        # remap regimes to new modelTypes
        mappingo <- attr(o, "mapping")
        for(r in seq_len(PCMNumRegimes(o))) {
          mappingo[r] <- m[mappingo[r]]
        }
        attr(o, "mapping") <- mappingo
        attr(o, "modelTypes") <- modelTypes
      }
    }
    o
  } else if(is.GaussianPCM(o)) {
    # We need to discover the corresponding MixedGaussian model type.
    pcmModelType <- class(o)[1L]
    typeParams <- strsplit(pcmModelType, split = "__", fixed = TRUE)[[1L]]
    modelType <- typeParams[1L]

    spec <- attr(o, "spec")
    paramNames <- names(spec)
    globalParams <- list()

    for(param in paramNames) {
      if(param == "X0" && !is.Global(spec[[param]])) {
        stop("as.MixedGaussian:: The parameter X0 should be global for all regimes.")
      }
      if(param == "Sigmae_x" && !is.Global(spec[[param]])) {
        globalParams[["Sigmae_x"]] <- structure(
          0.0,
          class = c("MatrixParameter", "_Omitted", "_Global"),
          description = "Global upper triangular factor of the non-phylogenetic variance-covariance")
      }
      if(is.Global(spec[[param]])) {
        globalParams[[param]] <- spec[[param]]
        modelType <- paste0(modelType, "__Omitted_", param)
      } else if(is.Omitted(spec[[param]])) {
        modelType <- paste0(modelType, "__Omitted_", param)
      } else {
        paramType <- do.call(paste0, as.list(class(spec[[param]])[-1L]))
        modelType <- paste0(modelType, "_", paramType, "_", param)
      }
    }

    regimes <- PCMRegimes(o)
    mapping <- rep(1L, length(regimes))
    names(mapping) <- as.character(regimes)

    if(!is.null(modelTypes)) {
      m <- match(modelType, modelTypes)
      if(is.na(m)) {
        stop(paste0("as.MixedGaussian:: ",
                    "The deduced Gaussian model type for the regimes (",
                    modelType, ") should match an element from modelTypes."))
      } else {
        mapping[] <- m
      }
    } else {
      modelTypes <- modelType
    }

    mg <- do.call(
      MixedGaussian,
      c(list(k = PCMNumTraits(o), modelTypes = modelTypes, mapping = mapping),
        globalParams))

    for(param in paramNames) {
      if(is.Global(spec[[param]]) && !is.Omitted(spec[[param]])) {
        mg[[param]][] <- o[[param]][]
      } else {
        for(regime in regimes) {
          if(is.ScalarParameter(spec[[param]])) {
            mg[[as.character(regime)]][[param]][] <- o[[param]][regime]
          } else if(is.VectorParameter(spec[[param]])) {
            mg[[as.character(regime)]][[param]][, 1L] <- o[[param]][, regime]
          } else if(is.MatrixParameter(spec[[param]])) {
            mg[[as.character(regime)]][[param]][,, 1L] <- o[[param]][,, regime]
          }
        }
      }
    }

    mg
  } else {
    stop(paste0(
      "as.MixedGaussian:: ",
      "This type of conversion is only possible if o is a GaussianPCM."))
  }
}

#' Create a template MixedGaussian object containing a regime for each model type
#'
#' @param mg a MixedGaussian object or an object that can be converted to such
#' via \code{\link{as.MixedGaussian}}.
#' @param modelTypes a (possibly named) character string
#' vector. Each such string denotes a mixed Gaussian regime model class, e.g.
#' the result of calling \code{MGPMDefaultModelTypes()}. If specified, an
#' attempt is made to match \code{PCMModelTypes(as.MixedGaussian(mg))}
#' with the elements of \code{modelTypes} and an error is raised if the match
#' fails. If not named, the model
#' types and regimes in the resulting MixedGaussian object are named by the
#' capital latin letters A,B,C,.... Default: \code{NULL}, which is interpreted
#' as \code{PCMModelTypes(as.MixedGaussian(mg, NULL))}.
#' @return a MixedGaussian with the same global parameter settings as for mg,
#' the same modelTypes as \code{modelTypes}, and with a regime for each model type.
#' The function will stop with an error if \code{mg} is not convertible to
#' a MixedGaussian object or if there is a mismatch between the model types in
#' \code{mg} and \code{modelTypes}.
#' @examples
#' mg <- MixedGaussianTemplate(PCMBaseTestObjects$model.ab.123.bSigmae_x)
#' mgTemplBMOU <- MixedGaussianTemplate(PCMBaseTestObjects$model.OU.BM)
#' @export
MixedGaussianTemplate <- function(mg, modelTypes = NULL) {
  mg <- as.MixedGaussian(mg, modelTypes = modelTypes)
  if(is.null(modelTypes)) {
    modelTypes <- PCMModelTypes(mg)
    if(is.null(names(modelTypes))) {
      names(modelTypes) <- LETTERS[seq_along(modelTypes)]
    }
  }
  k <- PCMNumTraits(mg)
  mapping <- structure(seq_along(modelTypes), names = names(modelTypes))

  spec <- attr(mg, "spec")
  globalParams <- spec[setdiff(names(spec), PCMRegimes(mg))]

  do.call(
    MixedGaussian,
    c(list(k = k, modelTypes = modelTypes, mapping = mapping), globalParams))
}

#' @export
PCMMapModelTypesToRegimes.MixedGaussian <- function(model, tree, ...) {
  uniqueRegimes <- unique(PCMRegimes(tree))
  res <- attr(model, "mapping", exact=TRUE)
  if(!is.null(names(res))) {
    res <- res[as.character(uniqueRegimes)]
  } else {
    names(res) <- as.character(uniqueRegimes)
  }
  res
}
