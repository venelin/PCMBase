# Copyright 2016-2019 Venelin Mitov, Krzysztof Bartoszek
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
PCMParentClasses.BM_drift <- function(model) {
  c("GaussianPCM", "PCM")
}

#' @export
PCMDescribe.BM_drift <- function(model, ...) {
  "Brownian motion model with drift"
}

#' @export
PCMCond.BM_drift <- function(
  tree, model, r = 1, metaI = PCMInfo(NULL, tree, model, verbose = verbose),
  verbose=FALSE) {


  Sigma_x <- if(is.Global(model$Sigma_x)){as.matrix(model$Sigma_x)} else {as.matrix(model$Sigma_x[,, r])}
  Sigma <- Sigma_x %*% t(Sigma_x)
  if(!is.null(model$Sigmae_x)) {
    Sigmae_x <- if(is.Global(model$Sigmae_x)){as.matrix(model$Sigmae_x)} else {as.matrix(model$Sigmae_x[,,r])}
    Sigmae <- Sigmae_x %*% t(Sigmae_x)
  } else {
    Sigmae <- NULL
  }

  if(!is.null(model$h_drift)) {
    h_drift <- if(is.Global(model$h_drift)) as.vector(model$h_drift) else model$h_drift[, r]
  }else{
    h_drift <- rep(0,nrow(Sigma_x))
  }

  V <- PCMCondVOU(matrix(0, nrow(Sigma), ncol(Sigma)), Sigma, Sigmae)
  omega <- function(t, edgeIndex, metaI) {
    t*h_drift
  }
  Phi <- function(t, edgeIndex, metaI, e_Ht = NULL) {
    diag(nrow(Sigma))
  }
  list(omega = omega, Phi = Phi, V = V)
}

#' @export
PCMDescribeParameters.BM_drift <- function(model, ...) {
  list(
    X0 = "trait values at the root",
    h_drift = "drift vector modyfing the expectation",
    Sigma_x = "Cholesky factor of the unit-time variance rate",
    Sigmae_x = "Cholesky factor of the non-heritable variance or the variance of the measurement error")
}

#' @export
PCMListParameterizations.BM_drift <- function(model, ...) {
  list(
    X0 = list(
      c("VectorParameter", "_Global"),
      c("VectorParameter", "_Fixed", "_Global"),
      c("VectorParameter", "_AllEqual", "_Global"),
      c("VectorParameter", "_Omitted")),
    h_drift = list(
      c("VectorParameter"),
      c("VectorParameter", "_Fixed"),
      c("VectorParameter", "_AllEqual"),
       c("VectorParameter", "_Omitted")),

    Sigma_x = list(
      c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal")),

    Sigmae_x = list(
      c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
      c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
      c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
      c("MatrixParameter", "_Omitted"))
  )
}

#' @export
PCMListDefaultParameterizations.BM_drift <- function(model, ...) {
  list(
    X0 = list(
      c("VectorParameter", "_Global"),
      c("VectorParameter", "_Omitted")
    ),
    h_drift = list(
      c("VectorParameter")),

    Sigma_x = list(
        c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
        c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
        c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal")
      ),

    Sigmae_x = list(
      c("MatrixParameter", "_Omitted"))
  )
}

#' @export
PCMSpecify.BM_drift <- function(model, ...) {
  spec <- list(
    X0 = structure(0.0, class = c('VectorParameter', '_Global'),
                   description = 'trait values at the root'),
    h_drift = structure(0.0, class = c('VectorParameter'),
                   description = 'drift vector modyfing the expectation'),
    Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
                        description = 'Cholesky factor of the unit-time variance rate'),
    Sigmae_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
                         description = 'Cholesky factor of the non-heritable variance or the variance of the measurement error'))
  attributes(spec) <- attributes(model)
  if(is.null(names(spec))) names(spec) <- c('X0', 'h_drift', 'Sigma_x', 'Sigmae_x')
  if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
  spec
}
