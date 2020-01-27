# This file was auto-generated through a call to PCMGenerateModelTypes()
# Do not edit by hand.


#' @export
PCMParentClasses.BM__Global_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("BM", "GaussianPCM", "PCM")

#' @export
PCMSpecify.BM__Global_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM__Omitted_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("BM", "GaussianPCM", "PCM")

#' @export
PCMSpecify.BM__Omitted_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM__Global_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("BM", "GaussianPCM", "PCM")

#' @export
PCMSpecify.BM__Global_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("BM", "GaussianPCM", "PCM")

#' @export
PCMSpecify.BM__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM__Global_X0__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("BM", "GaussianPCM", "PCM")

#' @export
PCMSpecify.BM__Global_X0__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM__Omitted_X0__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("BM", "GaussianPCM", "PCM")

#' @export
PCMSpecify.BM__Omitted_X0__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM__Global_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("BM", "GaussianPCM", "PCM")

#' @export
PCMSpecify.BM__Global_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM__Omitted_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("BM", "GaussianPCM", "PCM")

#' @export
PCMSpecify.BM__Omitted_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM__Global_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("BM", "GaussianPCM", "PCM")

#' @export
PCMSpecify.BM__Global_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("BM", "GaussianPCM", "PCM")

#' @export
PCMSpecify.BM__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM__Global_X0__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("BM", "GaussianPCM", "PCM")

#' @export
PCMSpecify.BM__Global_X0__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM__Omitted_X0__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("BM", "GaussianPCM", "PCM")

#' @export
PCMSpecify.BM__Omitted_X0__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM__Global_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("BM", "GaussianPCM", "PCM")

#' @export
PCMSpecify.BM__Global_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM__Omitted_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("BM", "GaussianPCM", "PCM")

#' @export
PCMSpecify.BM__Omitted_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM__Global_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("BM", "GaussianPCM", "PCM")

#' @export
PCMSpecify.BM__Global_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("BM", "GaussianPCM", "PCM")

#' @export
PCMSpecify.BM__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM__Global_X0__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("BM", "GaussianPCM", "PCM")

#' @export
PCMSpecify.BM__Global_X0__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.BM__Omitted_X0__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("BM", "GaussianPCM", "PCM")

#' @export
PCMSpecify.BM__Omitted_X0__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Global_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")

#' @export
PCMSpecify.OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root")),
H = structure(0.0, class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
description = c("adaptation rate matrix")),
Theta = structure(0.0, class = c("VectorParameter"),
description = c("long-term optimum")),
Sigma_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the unit-time variance rate")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
description = c("factor of the non-heritable variance or the variance of the measurement error")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.White__Global_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("White", "GaussianPCM", "PCM")

#' @export
PCMSpecify.White__Global_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root of the tree (for White model this is the mean vector)")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-phylogenetic variance-covariance matrix")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.White__Fixed_Global_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("White", "GaussianPCM", "PCM")

#' @export
PCMSpecify.White__Fixed_Global_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Fixed", "_Global"),
description = c("trait values at the root of the tree (for White model this is the mean vector)")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-phylogenetic variance-covariance matrix")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.White__AllEqual_Global_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("White", "GaussianPCM", "PCM")

#' @export
PCMSpecify.White__AllEqual_Global_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_AllEqual", "_Global"),
description = c("trait values at the root of the tree (for White model this is the mean vector)")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-phylogenetic variance-covariance matrix")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.White__Omitted_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("White", "GaussianPCM", "PCM")

#' @export
PCMSpecify.White__Omitted_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root of the tree (for White model this is the mean vector)")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-phylogenetic variance-covariance matrix")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.White__Global_X0__Diagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("White", "GaussianPCM", "PCM")

#' @export
PCMSpecify.White__Global_X0__Diagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root of the tree (for White model this is the mean vector)")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-phylogenetic variance-covariance matrix")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.White__Fixed_Global_X0__Diagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("White", "GaussianPCM", "PCM")

#' @export
PCMSpecify.White__Fixed_Global_X0__Diagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Fixed", "_Global"),
description = c("trait values at the root of the tree (for White model this is the mean vector)")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-phylogenetic variance-covariance matrix")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.White__AllEqual_Global_X0__Diagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("White", "GaussianPCM", "PCM")

#' @export
PCMSpecify.White__AllEqual_Global_X0__Diagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_AllEqual", "_Global"),
description = c("trait values at the root of the tree (for White model this is the mean vector)")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-phylogenetic variance-covariance matrix")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.White__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("White", "GaussianPCM", "PCM")

#' @export
PCMSpecify.White__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root of the tree (for White model this is the mean vector)")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-phylogenetic variance-covariance matrix")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.White__Global_X0__ScalarDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("White", "GaussianPCM", "PCM")

#' @export
PCMSpecify.White__Global_X0__ScalarDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root of the tree (for White model this is the mean vector)")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-phylogenetic variance-covariance matrix")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.White__Fixed_Global_X0__ScalarDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("White", "GaussianPCM", "PCM")

#' @export
PCMSpecify.White__Fixed_Global_X0__ScalarDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Fixed", "_Global"),
description = c("trait values at the root of the tree (for White model this is the mean vector)")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-phylogenetic variance-covariance matrix")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.White__AllEqual_Global_X0__ScalarDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("White", "GaussianPCM", "PCM")

#' @export
PCMSpecify.White__AllEqual_Global_X0__ScalarDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_AllEqual", "_Global"),
description = c("trait values at the root of the tree (for White model this is the mean vector)")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-phylogenetic variance-covariance matrix")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.White__Omitted_X0__ScalarDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model) c("White", "GaussianPCM", "PCM")

#' @export
PCMSpecify.White__Omitted_X0__ScalarDiagonal_WithNonNegativeDiagonal_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root of the tree (for White model this is the mean vector)")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
description = c("factor of the non-phylogenetic variance-covariance matrix")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.White__Global_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("White", "GaussianPCM", "PCM")

#' @export
PCMSpecify.White__Global_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root of the tree (for White model this is the mean vector)")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-phylogenetic variance-covariance matrix")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.White__Fixed_Global_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("White", "GaussianPCM", "PCM")

#' @export
PCMSpecify.White__Fixed_Global_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Fixed", "_Global"),
description = c("trait values at the root of the tree (for White model this is the mean vector)")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-phylogenetic variance-covariance matrix")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.White__AllEqual_Global_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("White", "GaussianPCM", "PCM")

#' @export
PCMSpecify.White__AllEqual_Global_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_AllEqual", "_Global"),
description = c("trait values at the root of the tree (for White model this is the mean vector)")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-phylogenetic variance-covariance matrix")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.White__Omitted_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("White", "GaussianPCM", "PCM")

#' @export
PCMSpecify.White__Omitted_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root of the tree (for White model this is the mean vector)")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-phylogenetic variance-covariance matrix")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.White__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("White", "GaussianPCM", "PCM")

#' @export
PCMSpecify.White__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root of the tree (for White model this is the mean vector)")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-phylogenetic variance-covariance matrix")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.White__Fixed_Global_X0__Diagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("White", "GaussianPCM", "PCM")

#' @export
PCMSpecify.White__Fixed_Global_X0__Diagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Fixed", "_Global"),
description = c("trait values at the root of the tree (for White model this is the mean vector)")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-phylogenetic variance-covariance matrix")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.White__AllEqual_Global_X0__Diagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("White", "GaussianPCM", "PCM")

#' @export
PCMSpecify.White__AllEqual_Global_X0__Diagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_AllEqual", "_Global"),
description = c("trait values at the root of the tree (for White model this is the mean vector)")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-phylogenetic variance-covariance matrix")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.White__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("White", "GaussianPCM", "PCM")

#' @export
PCMSpecify.White__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root of the tree (for White model this is the mean vector)")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-phylogenetic variance-covariance matrix")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.White__Global_X0__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("White", "GaussianPCM", "PCM")

#' @export
PCMSpecify.White__Global_X0__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Global"),
description = c("trait values at the root of the tree (for White model this is the mean vector)")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-phylogenetic variance-covariance matrix")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.White__Fixed_Global_X0__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("White", "GaussianPCM", "PCM")

#' @export
PCMSpecify.White__Fixed_Global_X0__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Fixed", "_Global"),
description = c("trait values at the root of the tree (for White model this is the mean vector)")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-phylogenetic variance-covariance matrix")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.White__AllEqual_Global_X0__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("White", "GaussianPCM", "PCM")

#' @export
PCMSpecify.White__AllEqual_Global_X0__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_AllEqual", "_Global"),
description = c("trait values at the root of the tree (for White model this is the mean vector)")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-phylogenetic variance-covariance matrix")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


#' @export
PCMParentClasses.White__Omitted_X0__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model) c("White", "GaussianPCM", "PCM")

#' @export
PCMSpecify.White__Omitted_X0__ScalarDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x <- function(model, ...) {
spec <- list(
X0 = structure(0.0, class = c("VectorParameter", "_Omitted"),
description = c("trait values at the root of the tree (for White model this is the mean vector)")),
Sigmae_x = structure(0.0, class = c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
description = c("factor of the non-phylogenetic variance-covariance matrix")))
attributes(spec) <- attributes(model)
if(is.null(names(spec))) names(spec) <- c("X0", "Sigmae_x")
if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
spec
}


