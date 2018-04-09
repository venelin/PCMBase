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

#' @name White
#' @title White Gaussian PCM ignoring phylogenetic history
#' @description White model ignoring phylogenetic history, treating trait values
#'  as independent samples from a k-variate Gaussian.
#' @details Calculating likelihoods for this model does not work if the global
#' option PCMBase.Singular.Skip is set to FALSE.
NULL

#' @inherit PCMParentClasses
#' @export
PCMParentClasses.White <- function(model) {
  c("GaussianPCM", "PCM")
}

#' @inherit PCMDescribe
#'
#' @export
PCMDescribe.White <- function(model, ...) {
  "White model ignoring phylogenetic history, treating trait values as independent samples
  from a k-variate Gaussian."
}

PCMInfo.White <- function(X, tree, model, verbose = FALSE) {
  res <- NextMethod()
  res$PCMBase.Skip.Singular <- TRUE
  res$PCMBase.Threshold.Skip.Singular <- Inf
  res
}

#' @export
PCMCond.White <- function(tree, model, r=1, metaI = PCMInfo(NULL, tree, model, verbose), verbose=FALSE) {
  Sigmae <- as.matrix(model$Sigmae_x[,,r]) %*% t(as.matrix(model$Sigmae_x[,,r]))

  V <- function(t, edgeIndex, metaI, e_Ht = NULL) {
    if(metaI$edge[edgeIndex,2] <= metaI$N) {
      Sigmae
    } else {
      matrix(as.double(0), metaI$k, metaI$k)
    }
  }

  omega <- function(t, edgeIndex, metaI) {
    rep(0, PCMNumTraits(model))
  }
  Phi <- function(t, edgeIndex, metaI, e_Ht = NULL) {
    diag(PCMNumTraits(model))
  }
  list(omega = omega, Phi = Phi, V = V)
}

#' @export
PCMSpecifyParams.White <- function(model, ...) {
  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  list(
    X0 = list(default = rep(0, k),
              type = c("gvector", "full"),
              description = "trait vector at the root; global for all model regimes"),
    Sigmae_x = list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "upper.tri.diag", "positive.diag"),
                  description = "variance-covariance matrix for the non-phylogenetic trait component"))
}

#' @export
PCMDescribe.White__posdiagSigmae_x <- function(model, ...) "White with positive diagonal Sigmae_x (i.e. uncorrelated traits)."
#' @export
PCMParentClasses.White__posdiagSigmae_x <- function(model) c("White", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.White__posdiagSigmae_x <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$Sigmae_x = list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "variance-covariance matrix for the non-phylogenetic trait component")
  spec
}

#' @export
PCMDescribe.White__noX0 <- function(model, ...) "White without X0."
#' @export
PCMParentClasses.White__noX0 <- function(model) c("White", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.White__noX0 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec[!sapply(spec, is.null)]
}


#' @export
PCMDescribe.White__noX0__posdiagSigmae_x <- function(model, ...) "White without X0 and with positive diagonal Sigmae_x (i.e. uncorrelated traits)."
#' @export
PCMParentClasses.White__noX0__posdiagSigmae_x <- function(model) c("White", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.White__noX0__posdiagSigmae_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)
  spec$Sigmae_x = list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "variance-covariance matrix for the non-phylogenetic trait component")
  spec[!sapply(spec, is.null)]
}


