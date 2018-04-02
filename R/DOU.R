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

#' @export
PCMParentClasses.DOU <- function(model) {
  c("GaussianPCM", "PCM")
}

#' @export
PCMCond.DOU <- function(tree, model, r=1, metaI = PCMInfo(NULL, tree, model, verbose), verbose=FALSE) {
  H1 <- as.matrix(model$H1[,,r])
  H2 <- as.matrix(model$H2[,,r])
  Theta <- model$Theta[,r]
  Sigma <- as.matrix(model$Sigma_x[,,r]) %*% t(as.matrix(model$Sigma_x[,,r]))
  if(!is.null(model$Sigmae)) {
    Sigmae <- as.matrix(model$Sigmae_x[,,r]) %*% t(as.matrix(model$Sigmae_x[,,r]))
  } else {
    Sigmae <- NULL
  }

  V <- PCMCondVOU(H2, Sigma, Sigmae, threshold.Lambda_ij = metaI$PCMBase.Threshold.Lambda_ij)
  omega <- function(t, edgeIndex, metaI, e_H1t = NULL) {
    if(is.null(e_H1t)) {
      e_H1t <- expm(-t*H1)
    }
    I <- diag(nrow(H1))
    (I-e_H1t) %*% Theta
  }
  Phi <- function(t, edgeIndex, metaI, e_H1t = NULL) {
    if(is.null(e_H1t)) {
      expm(-t*H1)
    } else {
      e_H1t
    }
  }
  list(omega = omega, Phi = Phi, V = V)
}

#' @export
PCMDescribe.DOU <- function(model, ...) {
  "Two-speed Ornstein-Uhlenbeck branching stochastic process model with a non-phylogenetic variance component"
}

#' @export
PCMSpecifyParams.DOU <- function(model, ...) {
  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  list(
    X0 = list(default = rep(0, k),
              type = c("gvector", "full"),
              description = "trait vector at the root; global for all model regimes"),
    H1 = list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
              type = c("matrix", "full"),
              description = "rate of convergence towards Theta"),
    H2 = list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
              type = c("matrix", "full"),
              description = "rate of decorrelation between pairs of tips"),
    Theta = list(default = array(0, dim = c(k, R), dimnames = list(NULL, regimes)),
                 type = c("vector", "full"),
                 description = "long-term optimum"),
    Sigma_x = list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                 type = c("matrix", "upper.tri.diag", "positive.diag"),
                 description = "Upper triangular Choleski factor of the unit-time variance-covariance matrix of the BM process"),
    Sigmae_x = list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "upper.tri.diag", "positive.diag"),
                  description = "Upper triangular Choleski factor of the variance-covariance matrix for the non-phylogenetic trait component"))
}



#' @export
PCMDescribe.DOU1 <- function(model, ...) "DOU without X0."
#' @export
PCMParentClasses.DOU1 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU1 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU2 <- function(model, ...) "DOU without Sigmae_x."
#' @export
PCMParentClasses.DOU2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU2 <- function(model, ...) {
  spec <- NextMethod()
  spec$Sigmae_x <- NULL
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU3 <- function(model, ...) "DOU without X0 and Sigmae_x."
#' @export
PCMParentClasses.DOU3 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU3 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL
  spec[!sapply(spec, is.null)]
}

