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
PCMParentClasses.BM <- function(model) {
  c("GaussianPCM", "PCM")
}

#' @export
PCMDescribe.BM <- function(model, ...) {
  "Brownian motion model"
}

#' @export
PCMCond.BM <- function(tree, model, r=1, metaI = PCMInfo(NULL, tree, model, verbose), verbose=FALSE) {

  Sigma_x <- if(is.Global(model$Sigma_x)) as.matrix(model$Sigma_x) else as.matrix(model$Sigma_x[,, r])
  Sigma <- Sigma_x %*% t(Sigma_x)
  if(!is.null(model$Sigmae_x)) {
    Sigmae_x <- if(is.Global(model$Sigmae_x)) as.matrix(model$Sigmae_x) else as.matrix(model$Sigmae_x[,,r])
    Sigmae <- Sigmae_x %*% t(Sigmae_x)
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

#' @export
PCMDescribeParameters.BM <- function(model, ...) {
  list(
    X0 = "trait values at the root",
    Sigma_x = "Choleski factor of the unit-time variance rate",
    Sigmae_x = "Choleski factor of the non-heritable variance or the variance of the measurement error")
}

#' @export
PCMListParameterizations.BM <- function(model, ...) {
  list(
    X0 = list(c("VectorParameter", "_Global"),
              c("VectorParameter", "_Fixed", "_Global"),
              c("VectorParameter", "_AllEqual", "_Global"),
              c("VectorParameter", "_Omitted")),
    Sigma_x = list(c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
                   c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
                   c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal")),

    Sigmae_x = list(c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
                    c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
                    c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
                    c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
                    c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
                    c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
                    c("MatrixParameter", "_Omitted"))
  )
}


#' @export
PCMSpecify.BM <- function(model, ...) {
  spec <- list(
    X0 = structure(0.0, class = c('VectorParameter', '_Global'),
                   description = 'trait values at the root'),
    Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
                        description = 'Choleski factor of the unit-time variance rate'),
    Sigmae_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
                         description = 'Choleski factor of the non-heritable variance or the variance of the measurement error'))
  attributes(spec) <- attributes(model)
  if(is.null(names(spec))) names(spec) <- c('X0', 'Sigma_x', 'Sigmae_x')
  if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
  spec
}
