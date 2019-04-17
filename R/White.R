# Copyright 2016-2019 Venelin Mitov
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

#' @export
PCMParentClasses.White <- function(model) {
  c("GaussianPCM", "PCM")
}

#' @export
PCMDescribe.White <- function(model, ...) {
  "White model"
}

#' @export
PCMInfo.White <- function(
  X, tree, model,
  SE = matrix(0.0, PCMNumTraits(model), PCMTreeNumTips(tree)),
  verbose = FALSE, preorder = NULL, ...) {
  if(is.Transformable(model)) {
    model <- PCMApplyTransformation(model)
  }

  res <- NextMethod()
  res$PCMBase.Skip.Singular <- TRUE
  res$PCMBase.Threshold.Skip.Singular <- Inf
  res
}

#' @export
PCMCond.White <- function(
  tree, model, r=1,
  metaI = PCMInfo(NULL, tree, model, verbose = verbose),
  verbose=FALSE) {
  if(!is.null(model$Sigmae_x)) {
    Sigmae_x <- if(is.Global(model$Sigmae_x)) as.matrix(model$Sigmae_x) else as.matrix(model$Sigmae_x[,,r])
    Sigmae <- Sigmae_x %*% t(Sigmae_x)
  } else {
    Sigmae <- NULL
  }

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
PCMDescribeParameters.White <- function(model, ...) {
  list(
    X0 = "trait values at the root of the tree (for White model this is the mean vector)",
    Sigmae_x = "Choleski factor of the non-phylogenetic variance-covariance matrix")
}

#' @export
PCMListParameterizations.White <- function(model, ...) {
  list(
    X0 = list(c("VectorParameter", "_Global"),
              c("VectorParameter", "_Fixed", "_Global"),
              c("VectorParameter", "_AllEqual", "_Global"),
              c("VectorParameter", "_Omitted")),
    Sigmae_x = list(c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
                    c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
                    c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
                    c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
                    c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
                    c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"))
  )
}

#' @export
PCMSpecify.White <- function(model, ...) {
  spec <- list(
    X0 = structure(0.0, class = c('VectorParameter', '_Global'),
                   description = 'trait values at the root of the tree (for White model this is the mean vector)'),
    Sigmae_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
                         description = 'Choleski factor of the non-phylogenetic variance-covariance matrix'))
  attributes(spec) <- attributes(model)
  if(is.null(names(spec))) names(spec) <- c('X0', 'Sigmae_x')
  if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
  spec
}
