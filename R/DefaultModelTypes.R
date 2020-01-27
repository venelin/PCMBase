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

#' Class names for the the default PCM and MGPM model types
#' @description Utility functions returning named character vector of the
#' model class-names for the default model types used for PCM and MixedGaussian
#' model construction.
#' @return Both, \code{PCMFDefaultModelTypes} and
#' \code{MGPMDefaultModelTypes} return a character string vector with
#' named elements (A,B,C,D,E,F) defined as follows
#' (Mitov et al. 2019a):
#' \describe{
#' \item{A. }{BM (H = 0, diagonal \eqn{\Sigma}): BM, uncorrelated traits.}
#' \item{B. }{BM (H = 0, symmetric \eqn{\Sigma}): BM, correlated traits.}
#' \item{C. }{OU (diagonal H, diagonal \eqn{\Sigma}): OU, uncorrelated traits.}
#' \item{D. }{OU (diagonal H, symmetric \eqn{\Sigma}): OU, correlated traits, but simple
#' (diagonal) selection strength matrix.}
#' \item{E. }{OU (symmetric H, symmetric \eqn{\Sigma}): An OU with nondiagonal symmetric H
#' and nondiagonal symmetric \eqn{\Sigma}.}
#' \item{F. }{OU (asymmetric H, symmetric \eqn{\Sigma}): An OU with nondiagonal asymmetric
#' H and nondiagonal symmetric \eqn{\Sigma}.}
#' }
#' The only difference between the two functions is that the model
#' types returned by \code{PCMFDefaultModelTypes} have a global
#' parameter X0, while the model types returned by
#' \code{MGPMFDefaultModelTypes} have an omitted parameter X0.
#'
#' @seealso Args_MixedGaussian_MGPMDefaultModelTypes
#' @references
#' [Mitov et al. 2019a] Mitov, V., Bartoszek, K., & Stadler, T. (2019). Automatic generation of
#'  evolutionary hypotheses using mixed Gaussian phylogenetic models.
#'  Proceedings of the National Academy of Sciences of the United States of
#'  America, 35, 201813823. http://doi.org/10.1073/pnas.1813823116
#'
#' @export
PCMDefaultModelTypes <- function() {
  c(
    # BM; independent traits
    A = "BM__Global_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
    # BM; dependent traits
    B = "BM__Global_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
    # OU; independent traits (diagonal H and diagonal Sigma)
    C = "OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
    # OU; dependtent traits (diagonal H and non-diagonal Sigma)
    D = "OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
    # OU; dependent traits (Symmetric H and non-diagonal Sigma)
    E = "OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
    # OU; dependent traits with causality (non-symmetric H and non-diagonal Sigma)
    F = "OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"
  )
}

#' @rdname PCMDefaultModelTypes
#' @export
MGPMDefaultModelTypes <- function() {
  c(
    # BM; independent traits
    A = "BM__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
    # BM; dependent traits
    B = "BM__Omitted_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
    # OU; independent traits (diagonal H and diagonal Sigma)
    C = "OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
    # OU; dependtent traits (diagonal H and non-diagonal Sigma)
    D = "OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
    # OU; dependent traits (Symmetric H and non-diagonal Sigma)
    E = "OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
    # OU; dependent traits with causality (non-symmetric H and non-diagonal Sigma)
    F = "OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"
  )
}

#' Arguments to be passed to the constructor MixedGaussian when constructing
#' a MGPM model with some of the default MGPM model types.
#' @param omitGlobalSigmae_x logical, indicating if the returned list should specify
#' the global Sigmae_x parameter as '_Omitted'. Default: TRUE.
#' @seealso MGPMDefaultModelTypes
#' @return a list of named arguments. Currently only a named element Sigmae_x
#' with specification depending on \code{omitGlobalSigmae_x}.
#' @export
Args_MixedGaussian_MGPMDefaultModelTypes <- function(omitGlobalSigmae_x = TRUE) {
  if(omitGlobalSigmae_x) {
    list(
      Sigmae_x = structure(
        0.0, class = c("MatrixParameter", "_Omitted"),
        description = "upper triangular factor of the non-phylogenetic variance-covariance")
    )
  } else {
    list(Sigmae_x = structure(
      0.0,
      class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
      description = "Upper triangular factor of the non-phylogenetic variance-covariance")
    )
  }
}





#' Class name for the scalar OU MGPM model type
#' @return a character vector of one named element (ScalarOU)
#' @export
MGPMScalarOUType <- function() {
  c(
    ScalarOU = "OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"
  )
}

#' Arguments for the MixedGaussian constructor for scalar OU MGPM models.
#' @return a list.
#' @export
Args_MixedGaussian_MGPMScalarOUType <- function() {
  list(
    H = structure(
      0.0,
      class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
      description = "adaptation rate matrix"),
    Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
                         description = "upper triangular factor of the non-phylogenetic variance-covariance")
  )
}

#' Class name for the SURFACE OU MGPM model type
#' @return a character vector of one named element (SURFACE)
#' @export
MGPMSurfaceOUType <- function() {
  c(
    SURFACE = "OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x"
  )
}



#' Arguments for the MixedGaussian constructor for SURFACE OU MGPM models.
#' @return a list.
#' @export
Args_MixedGaussian_MGPMSurfaceOUType <- function() {
  list(
    H = structure(
      0.0,
      class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
      description = "adaptation rate matrix"),
    Sigma_x = structure(
      0.0,
      class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
      description = "unit-time variance parameter of the OU-process"),
    Sigmae_x = structure(
      0.0, class = c("MatrixParameter", "_Omitted"),
      description = "upper triangular factor of the non-phylogenetic variance-covariance")
  )
}

# The folloaing line of code is executed manually, while the current working
# directory is <package-home>/R/.
# PCMGenerateModelTypes(source = "AutoGeneratedModelTypes.R")
# This will generate the R-file
# AutoGeneratedModelTypes.R. Once this is done Rerun Roxygen2 to update the
# NAMESPACE file; rebuild the package.

