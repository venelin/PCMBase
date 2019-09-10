# This file was auto-generated through a call to PCMGenerateModelTypes()
# Do not edit by hand.


#' @export
PCMParentClasses.BM__Global_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('BM', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.BM__Global_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM__Omitted_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('BM', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.BM__Omitted_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM__Global_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('BM', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.BM__Global_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('BM', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.BM__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM__Global_X0__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('BM', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.BM__Global_X0__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM__Omitted_X0__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('BM', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.BM__Omitted_X0__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM_drift__Global_X0__h_drift__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('BM_drift', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.BM_drift__Global_X0__h_drift__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
h_drift = structure(0.0, class = c('VectorParameter'),
description = 'drift vector modyfing the expectation'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
description = 'Upper triangular factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'Upper triangular factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'h_drift', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM_drift__Omitted_X0__h_drift__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('BM_drift', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.BM_drift__Omitted_X0__h_drift__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
h_drift = structure(0.0, class = c('VectorParameter'),
description = 'drift vector modyfing the expectation'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
description = 'Upper triangular factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'Upper triangular factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'h_drift', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM_drift__Global_X0__h_drift__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('BM_drift', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.BM_drift__Global_X0__h_drift__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
h_drift = structure(0.0, class = c('VectorParameter'),
description = 'drift vector modyfing the expectation'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'Upper triangular factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'Upper triangular factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'h_drift', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM_drift__Omitted_X0__h_drift__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('BM_drift', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.BM_drift__Omitted_X0__h_drift__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
h_drift = structure(0.0, class = c('VectorParameter'),
description = 'drift vector modyfing the expectation'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'Upper triangular factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'Upper triangular factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'h_drift', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM_drift__Global_X0__h_drift__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('BM_drift', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.BM_drift__Global_X0__h_drift__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
h_drift = structure(0.0, class = c('VectorParameter'),
description = 'drift vector modyfing the expectation'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal'),
description = 'Upper triangular factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'Upper triangular factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'h_drift', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM_drift__Omitted_X0__h_drift__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('BM_drift', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.BM_drift__Omitted_X0__h_drift__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
h_drift = structure(0.0, class = c('VectorParameter'),
description = 'drift vector modyfing the expectation'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal'),
description = 'Upper triangular factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'Upper triangular factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'h_drift', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_Diagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_Diagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_Diagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_Diagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_Diagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_Diagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_Diagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_Diagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_Diagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_Diagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_Diagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_Diagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_Diagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_Diagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_Diagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_Diagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_Diagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_Diagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_Diagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_Diagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_Diagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_Diagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Transformable'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Diagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_Diagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_Diagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Global'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c('OU', 'GaussianPCM', 'PCM')

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c('VectorParameter', '_Omitted'),
description = 'trait values at the root'),
H = structure(0.0, class = c('MatrixParameter', '_Schur', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Transformable', '_Global'),
description = 'adaptation rate matrix'),
Theta = structure(0.0, class = c('VectorParameter'),
description = 'long-term optimum'),
Sigma_x = structure(0.0, class = c('MatrixParameter', '_ScalarDiagonal', '_WithNonNegativeDiagonal', '_Global'),
description = 'factor of the unit-time variance rate'),
Sigmae_x = structure(0.0, class = c('MatrixParameter', '_Omitted'),
description = 'factor of the non-heritable variance or the variance of the measurement error'))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c('X0', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


