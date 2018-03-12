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

#' @name BM
#' @title Brownian motion PCMs
NULL

#' @inherit PCMParentClasses
#' @export
PCMParentClasses.BM <- function(model) {
  c("GaussianPCM", "PCM")
}

#' @inherit PCMDescribe
#'
#' @export
PCMDescribe.BM <- function(model, ...) {
  "Brownian motion (BM) branching stochastic process model with a non-phylogenetic variance component"
}

#' @export
PCMCond.BM <- function(tree, model, r=1, metaI = PCMInfo(NULL, tree, model, verbose), verbose=FALSE) {
  Sigma <- as.matrix(model$Sigma[,,r])
  if(!is.null(model$Sigmae)) {
    Sigmae <- as.matrix(model$Sigmae[,,r])
  } else {
    Sigmae <- NULL
  }

  V <- PCMCondVOU(matrix(0, nrow(Sigma), ncol(Sigma)), Sigma, Sigmae)
  omega <- function(t, edgeIndex, metaI) {
    rep(0, nrow(Sigma))
  }
  Phi <- function(t, edgeIndex, metaI, e_Ht = NULL) {
    diag(nrow(Sigma))
  }
  list(omega = omega, Phi = Phi, V = V)
}


#' @inherit PCMSpecifyParams
#' @export
PCMSpecifyParams.BM <- function(model, ...) {
  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  list(
    X0 = list(default = rep(0, k),
              type = c("gvector", "full"),
              description = "trait vector at the root; global for all model regimes"),
    Sigma = list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                 type = c("matrix", "symmetric"),
                 description = "unit-time variance-covariance matrix of the BM-process"),
    Sigmae = list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "symmetric"),
                  description = "variance-covariance matrix for the non-phylogenetic trait component"))
}

#' @export
PCMDescribe.BM1 <- function(model, ...) "BM without X0."
#' @export
PCMParentClasses.BM1 <- function(model) c("BM", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.BM1 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec[!sapply(spec, is.null)]
}
