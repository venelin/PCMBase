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

#' Conditional multivariate OU distribution
#' @inheritParams PCMCond
#'
#' @return a named list as specified in \code{\link{PCMCond}}
#' @importFrom expm expm
#' @export
PCMCond.OU <- function(tree, model, r=1, metaI=PCMInfo(NULL, tree, model, verbose), verbose=FALSE) {
  H <- as.matrix(model$H[,, r])
  Theta <- model$Theta[, r]
  Sigma <- as.matrix(model$Sigma[,,r])

  V <- PCMCondVOU(H, Sigma)
  omega <- function(t, edgeIndex, e_Ht = NULL) {
    if(is.null(e_Ht)) {
      e_Ht <- expm(-t*H)
    }
    I <- diag(nrow(H))
    (I-e_Ht) %*% Theta
  }
  Phi <- function(t, edgeIndex, e_Ht = NULL) {
    if(is.null(e_Ht)) {
      expm(-t*H)
    } else {
      e_Ht
    }
  }
  list(omega = omega, Phi = Phi, V = V)
}

#' @describeIn PCMDescribe
#' @export
PCMDescribe.OU <- function(model, ...) {
  "Ornstein-Uhlenbeck branching stochastic process model and a non-phylogenetic variance component"
}

#' @describeIn PCMSpecifyParams
#' @export
PCMSpecifyParams.OU <- function(model, ...) {
  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  list(
    X0 = list(default = rep(0, k),
              type = c("gvector", "full"),
              description = "trait vector at the root; global for all model regimes"),
    H = list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
             type = c("matrix", "full"),
             description = "Selection strength matirx"),
    Theta = list(default = array(0, dim = c(k, R), dimnames = list(NULL, regimes)),
                 type = c("vector", "full"),
                 description = "long-term optimum trait values"),
    Sigma = list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                 type = c("matrix", "symmetric"),
                 description = "unit-time variance-covariance matrix of the BM-process"),
    Sigmae = list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "symmetric"),
                  description = "variance-covariance matrix for the non-phylogenetic trait component"))
}
