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

#' @name MixedGaussian
#' @title Multiple regime Gaussian PCMs

#' @export
is.MixedGaussian <- function(x) inherits(x, "MixedGaussian")

#' @export
PCMParentClasses.MixedGaussian <- function(model) {
  c("GaussianPCM", "PCM")
}

#' @export
PCMDescribe.MixedGaussian <- function(model, ...) {
  "Mixed Gaussian model"
}

#' @export
PCMCond.MixedGaussian <- function(tree, model, r=1, metaI = PCMInfo(NULL, tree, model, verbose), verbose=FALSE) {
  if(! is.MixedGaussian(model) ) {
    stop("ERR:02501:PCMBase:MixedGaussian.R:PCMCond.MixedGaussian:: model should inherit from S3 class 'MixedGaussian'.")
  }
  if(!is.character(r)) {
    # integer regime number must be mapped to a character regime because
    # the model object may contain global parameter vectors and matrices apart
    # from models associated with regimes.
    r <- as.character(attr(model, "regimes", exact=TRUE)[r])
  }

  if(!is.null(model$Sigmae_x)) {
    Sigmae <- as.matrix(model$Sigmae_x) %*% t(as.matrix(model$Sigmae_x))
  } else {
    Sigmae <- NULL
  }

  if(is.null(Sigmae)) {
    # Sigmae_x is omitted in the parent model, so just return the result from the
    # sub-model
    PCMCond(tree, model[[r]], 1, metaI, verbose)
  } else {
    OmegaPhiV2 <- OmegaPhiV <- PCMCond(tree, model[[r]], 1, metaI, verbose)
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
PCMParamCount.MixedGaussian <- function(o, countRegimeChanges = FALSE, countModelTypes = FALSE,  offset = 0, k = 1, R = 1, parentModel = NULL) {
  p <- NextMethod()
  if(countModelTypes) {
    if(length(attr(o, "modelTypes", exact = TRUE)) > 1) {
      # + 1 has already been counted by the default PCM implementatioin (call
      # to NextMethod() above)
      p <- p + PCMNumRegimes(o) - 1
    }
  }
  p
}


#' Create a multi-regime Gaussian model (MixedGaussian)
#' @param k integer defining the number of traits.
#' @param modelTypes a character string vector with the class names of the
#' model-types that can possibly be included (assigned to regimes) in the MixedGaussian,
#' e.g. c("BM3", "OU3").
#' @param mapping a character string vector with elements from modelTypes or an
#' integer vector with elements between 1 and length(modelTypes)
#' mapping modelTypes to regimes, e.g. if \code{modelTypes = c("BM3", "OU3")} and
#' \code{mapping = c(a = 1, b = 1, c = 2, d = 1)} defines an MixedGaussian with four
#' different regimes with model-types BM3, BM3, OU3 and BM3, corresponding to each
#' regime. \code{mapping} does not have to be a named vector. If it is a named
#' vector, then all the names must correspond to valid regime names in a tree to
#' which the model will be fit (member tree$edge.regime should be a character vector).
#' If it is not a named vector then the positions of the elements correspond to the
#' regimes in their order given by the function \code{\link{PCMTreeUniqueRegimes}}
#' called on a tree object.
#' @param className a character string definingn a valid S3 class name for the
#' resulting MixedGaussian object. If not specified, a className is generated using the
#' expression \code{ paste0("MixedGaussian_", do.call(paste0, as.list(mapping)))}.
#' @param X0 specification for the global vector X0 to be used by all
#' models in the MixedGaussian.
#' @param ... specifications for other _Global parameters coming after X0.
#' @param Sigmae_x sepcification of a _Global Sigmae_x parameter. This is used
#' by Submodels only if they have Sigmae_x _Omitted.
#' @return an object of S3 class className inheriting from MixedGaussian, GaussianPCM and
#' PCM.
#'
#' @details If X0 is not NULL it has no sense to use model-types including X0 as a
#' parameter (e.g. use BM1 or BM3 insted of BM or BM2). Similarly if Sigmae_x is
#' not NULL there is no meaning in using model-types including Sigmae_x as a parameter,
#' (e.g. use OU2 or OU3 instead of OU or OU1).
#' @seealso \code{\link{PCMTreeUniqueRegimes}}
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
  Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
                       description = "Upper triangular Choleski factor of the non-phylogenetic variance-covariance")) {

  regimes <- if(is.null(names(mapping))) seq_len(length(mapping)) else names(mapping)

  if(is.character(mapping)) {
    mapping2 <- match(mapping, modelTypes)
    if(any(is.na(mapping2))) {
      stop(paste0("ERR:02511:PCMBase:MixedGaussian.R:MixedGaussian:: some of the models in mapping not found in modelTypes: ",
                  "modelTypes = ", toString(modelTypes), ", mapping =", toString(mapping)))
    } else {
      mapping <- mapping2
    }
  }

  mappingModelRegime <- modelTypes[mapping]

  spec <- list(X0 = X0, ...)

  for(m in 1:length(mapping)) {
    subModel <- structure(0.0, class = mappingModelRegime[m])
    class(subModel) <- c(class(subModel), PCMParentClasses(subModel))
    spec[[as.character(regimes[m])]] <- PCMSpecify(subModel)
    attr(spec[[as.character(regimes[m])]], "k") <- k
    attr(spec[[as.character(regimes[m])]], "regimes") <- 1
  }

  spec[["Sigmae_x"]] <- Sigmae_x

  #spec[] <- spec[!sapply(spec, is.null)]

  attr(spec, "k") <- k

  class(spec) <- c(className, "MixedGaussian", "GaussianPCM", "PCM")

  if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')

  res <- PCM(model = class(spec), modelTypes, k = k, regimes = regimes, spec = spec)
  attr(res, "mapping") <- mapping
  attr(res, "modelTypes") <- modelTypes
  attr(res, "spec") <- spec
  res
}


#' @export
PCMMapModelTypesToRegimes.MixedGaussian <- function(model, tree, ...) {
  uniqueRegimes <- PCMTreeUniqueRegimes(tree)
  res <- attr(model, "mapping", exact=TRUE)
  if(!is.null(names(res))) {
    res <- res[as.character(uniqueRegimes)]
  } else {
    names(res) <- as.character(uniqueRegimes)
  }
  res
}