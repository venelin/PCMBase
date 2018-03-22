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

#' @name MRG
#' @title Multiple regime Gaussian PCMs

#' @export
is.MRG <- function(x) inherits(x, "MRG")

#' @export
PCMParentClasses.MRG <- function(model) {
  c("GaussianPCM", "PCM")
}

#' @export
PCMDescribe.MRG <- function(model, ...) {
  "Multiple regime Gaussian model"
}

#' @export
PCMCond.MRG <- function(tree, model, r=1, metaI = PCMInfo(NULL, tree, model, verbose), verbose=FALSE) {
  if(! is.MRG(model) ) {
    stop("ERR:02501:PCMBase:MRG.R:PCMCond.MRG:: model should inherit from S3 class 'MRG'.")
  }
  if(!is.character(r)) {
    # integer regime number must be mapped to a character regime because
    # the model object may contain global parameter vectors and matrices apart
    # from models associated with regimes.
    r <- attr(model, "regimes")[r]
  }
  PCMCond(tree, model[[r]], 1, metaI, verbose)
}

#' @export
PCMGetVecParamsFull.MRG <- function(model, ...) {
  specParams <- attr(model, "specParams", exact = TRUE)
  R <- PCMNumRegimes(model)

  types <- sapply(names(model), function(name) specParams[[name]]$type[1])

  # the order of types is important!
  # gscalars, gvectors, gmatrices, which appear before the first model are
  # prepended to each model; gscalars, gvectors, gmatrices, which appear after the last model
  # are appended to each model.
  # the vector gmask below is TRUE for prependable g-params and FALSE for models and append-able g-params
  gmask <- !(cumsum(as.integer(types=="model")))
  names(gmask) <- names(model)

  gprep <- c()
  gapp <- c()
  for(name in names(model)) {
    if(gmask[name]) {
      gprep <- c(gprep, as.vector(model[[name]]))
    } else if(types[[name]] != "model") {
      gapp <- c(gapp, as.vector(model[[name]]))
    }
  }

  res <- do.call(c, lapply(names(model), function(name) {
    if(specParams[[name]]$type[1]=="model") {
      c(gprep, PCMGetVecParamsFull(model[[name]], ...), gapp)
    }
  }))
  unname(res)
}


#' Create a multi-regime Gaussian model (MRG)
#' @param k integer defining the number of traits.
#' @param models a character string vector with the class names of the possible
#' models included in the MRG, e.g. c("BM3", "OU3").
#' @param mapping a character string vector with elements from models or an
#' integer vector with elements between 1 and length(models)
#' mapping models to regimes (the positions in the mapping vector), e.g.
#' c(a = 1, b = 1, c = 2, d = 1) defines an MRG with four different regimes with
#' models BM3, BM3, OU3 and BM3, corresponding to each regime.
#' @param className a character string definingn a valid class name for the
#' resulting MRG object.
#' @param X0 NULL or a list defining a global vector X0 to be used by all
#' models in the MRG.
#' @param Sigmae NULL or a list defining a global matrix to be used as element
#' Sigmae by all models in the MRG.
#' @return an object of S3 class className inheriting from MRG, GaussianPCM and
#' PCM.
#'
#' @details If X0 is not NULL it has no sense to use models including X0 as a
#' parameter (e.g. use BM1 or BM3 insted of BM or BM2). Similarly if Sigmae is
#' not NULL there is no meaning in using models including Sigmae as a parameter,
#' (e.g. use OU2 or OU3 instead of OU or OU1).
#' @export
MRG <- function(
  k,
  models,
  mapping,
  className = paste0("MRG_", do.call(paste0, as.list(mapping))),
  X0 = list(default = rep(0, k),
            type = c("gvector", "full"),
            description = "trait vector at the root; global for all regimes"),
  Sigmae = list(default = matrix(0, nrow = k, ncol = k),
                type = c("gmatrix", "symmetric"),
                description = "variance-covariance matrix for the non-phylogenetic trait component")) {

  regimes <- if(is.null(names(mapping))) as.character(1:length(mapping)) else names(mapping)

  if(is.character(mapping)) {
    mapping2 <- match(mapping, models)
    if(any(is.na(mapping2))) {
      stop(paste0("ERR:02511:PCMBase:MRG.R:MRG:: some of the model-names in mapping not found in models: ",
                  "models = ", toString(models), ", mapping =", toString(mapping)))
    } else {
      mapping <- mapping2
    }
  }

  mappingModelRegime <- models[mapping]

  specParams <- list(X0 = X0)

  for(m in 1:length(mapping)) {
    specParams[[regimes[m]]] <- list(default = PCM(mappingModelRegime[m], k, 1), type = "model")
  }

  specParams[["Sigmae"]] <- Sigmae
  specParams <- specParams[!sapply(specParams, is.null)]

  PCM(model = c(className, "MRG", "GaussianPCM", "PCM"), k = k, regimes = regimes, specParams = specParams)
}

#' @export
PCMLowerBound.MRG <- function(model, X0 = NULL, Sigmae = NULL, ...) {
  model <- NextMethod()
  specParams <- attr(model, "specParams", exact = TRUE)

  if(!is.null(specParams$X0) && !is.null(X0)) {
    model$X0 <- X0
  }
  for(name in names(specParams)) {
    if(specParams[[name]]$type[1]=="model") {
      model[[name]] <- PCMLowerBound(model[[name]], X0 = X0, Sigmae = Sigmae, ...)
    }
  }
  if(!is.null(specParams$Sigmae) && !is.null(Sigmae)) {
    model$Sigmae <- Sigmae
  }
  model
}

#' @export
PCMUpperBound.MRG <- function(model, X0 = NULL, Sigmae = NULL, ...) {
  model <- NextMethod()
  specParams <- attr(model, "specParams", exact = TRUE)

  if(!is.null(specParams$X0) && !is.null(X0)) {
    model$X0 <- X0
  }
  for(name in names(specParams)) {
    if(specParams[[name]]$type[1]=="model") {
      model[[name]] <- PCMUpperBound(model[[name]], X0 = X0, Sigmae = Sigmae, ...)
    }
  }
  if(!is.null(specParams$Sigmae) && !is.null(Sigmae)) {
    model$Sigmae <- Sigmae
  }
  model
}

