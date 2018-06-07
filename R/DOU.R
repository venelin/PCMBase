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
  H1 <- if(is.Global(model$H1)) as.matrix(model$H1) else as.matrix(model$H1[,, r])

  H2 <- if(is.Global(model$H2)) as.matrix(model$H2) else as.matrix(model$H2[,, r])

  Theta <- if(is.Global(model$Theta)) as.vector(model$Theta) else model$Theta[, r]

  Sigma_x <- if(is.Global(model$Sigma_x)) as.matrix(model$Sigma_x) else as.matrix(model$Sigma_x[,, r])

  Sigma <- Sigma_x %*% t(Sigma_x)

  if(!is.null(model$Sigmae_x)) {
    Sigmae_x <- if(is.Global(model$Sigmae_x)) as.matrix(model$Sigmae_x) else as.matrix(model$Sigmae_x[,,r])
    Sigmae <- Sigmae_x %*% t(Sigmae_x)
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
  "Double-rate Ornstein-Uhlenbeck model"
}

#' @export
PCMDescribeParameters.DOU <- function(model, ...) {
  list(
    X0 = "trait values at the root",
    H1 = "adaptation rate matrix",
    H2 = "decorrelation rate matrix",
    Theta = "long-term optimum",
    Sigma_x = "Choleski factor of the unit-time variance rate",
    Sigmae_x = "Choleski factor of the non-heritable variance or the variance of the measurement error")
}

#' @export
PCMListParameterizations.DOU <- function(model, ...) {
  list(
    X0 = list(c("VectorParameter", "_Global"),
              c("VectorParameter", "_Fixed", "_Global"),
              c("VectorParameter", "_AllEqual", "_Global"),
              c("VectorParameter", "_Omitted")),
    H1 = list(c("MatrixParameter"),
              c("MatrixParameter", "_Fixed"),
              c("MatrixParameter", "_Symmetric"),
              c("MatrixParameter", "_Diagonal"),
              c("MatrixParameter", "_ScalarDiagonal"),
              c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
              c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
              c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
              c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
              c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
              c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
              c("MatrixParameter", "_Global"),
              c("MatrixParameter", "_Fixed", "_Global"),
              c("MatrixParameter", "_Symmetric", "_Global"),
              c("MatrixParameter", "_Diagonal", "_Global"),
              c("MatrixParameter", "_ScalarDiagonal", "_Global"),
              c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
              c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
              c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
              c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
              c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
              c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global")),
    H2 = list(c("MatrixParameter"),
              c("MatrixParameter", "_Fixed"),
              c("MatrixParameter", "_Symmetric"),
              c("MatrixParameter", "_Diagonal"),
              c("MatrixParameter", "_ScalarDiagonal"),
              c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
              c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
              c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
              c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
              c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
              c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
              c("MatrixParameter", "_Global"),
              c("MatrixParameter", "_Fixed", "_Global"),
              c("MatrixParameter", "_Symmetric", "_Global"),
              c("MatrixParameter", "_Diagonal", "_Global"),
              c("MatrixParameter", "_ScalarDiagonal", "_Global"),
              c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
              c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
              c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
              c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
              c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
              c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global")),

    Theta = list(c("VectorParameter"),
                 c("VectorParameter", "_Fixed"),
                 c("VectorParameter", "_AllEqual")),

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
PCMSpecify.DOU <- function(model, ...) {
  spec <- list(
    X0 = structure(0.0, class = c('VectorParameter', '_Global'),
                   description = 'trait values at the root'),
    H1 = structure(0.0, class = c('MatrixParameter'),
                   description = 'adaptation rate matrix'),
    H2 = structure(0.0, class = c('MatrixParameter'),
                   description = 'decorrelation rate matrix'),
    Theta = structure(0.0, class = c('VectorParameter'),
                      description = 'long-term optimum'),
    Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
                        description = 'Choleski factor of the unit-time variance rate'),
    Sigmae_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
                         description = 'Choleski factor of the non-heritable variance or the variance of the measurement error'))
  attributes(spec) <- attributes(model)
  if(is.null(names(spec))) names(spec) <- c('X0', 'H1', 'H2', 'Theta', 'Sigma_x', 'Sigmae_x')
  if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
  spec
}
