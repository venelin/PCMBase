#' Class names for the the default PCM and MGPM model types
#' @description Utility functions returning named character vector of the
#' model class-names for the default model types used for PCM and MixedGaussian
#' model construction.
#' @return a character string vector with named elements
#' (currently  A,B,C,D,E,F).
#' @seealso Args_MixedGaussian_MGPMDefaultModelTypes
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
#' @seealso MGPMDefaultModelTypes
#' @return a list of named arguments.
#' @export
Args_MixedGaussian_MGPMDefaultModelTypes <- function() {
  list(
    Sigmae_x = structure(
      0.0, class = c("MatrixParameter", "_Omitted"),
      description = "upper triangular Choleski factor of the non-phylogenetic variance-covariance")
  )
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
                         description = "upper triangular Choleski factor of the non-phylogenetic variance-covariance")
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
      description = "upper triangular Choleski factor of the non-phylogenetic variance-covariance")
  )
}

# GenerateDefaultParameterizations <- function(libname, pkgname){
#   PCMGenerateParameterizations(
#     structure(0.0, class = "OU"),
#     PCMTableParameterizations(structure(0.0, class = "OU"))[1])
#   PCMGenerateParameterizations(
#     structure(0.0, class = "DOU"),
#     PCMTableParameterizations(structure(0.0, class = "DOU"))[1])
#   PCMGenerateParameterizations(
#     structure(0.0, class = "BM"),
#     PCMTableParameterizations(structure(0.0, class = "BM"))[1])
#   PCMGenerateParameterizations(
#     structure(0.0, class = "JOU"),
#     PCMTableParameterizations(structure(0.0, class = "JOU"))[1])
#   PCMGenerateParameterizations(
#     structure(0.0, class = "White"),
#     PCMTableParameterizations(structure(0.0, class = "White"))[1])
# }
#