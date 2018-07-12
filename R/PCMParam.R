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


## API for loading or storing parameters from/to a vector ----------------------
# The PCM framework has a helper API for storing, transforming and interpreting
# model parameters. This is build upon a hierarchy of S3-classes. There are
# three types of S3 classes:
# 1. Storage classes:
#   - ScalarParameter: numeric parameter (a single number), local to a regime;
#   - VectorParameter: numeric vector of length k, local to a regime;
#   - MatrixParameter: numeric matrix of dimension k x k, local to a regime;
# 2.

##  scope
#' @export
is.Local <- function(o) { !inherits(o, "_Global") }
#' @export
is.Global <- function(o) { inherits(o, "_Global") }

## data-type
#' @export
is.ScalarParameter <- function(o) { inherits(o, "ScalarParameter") }
#' @export
is.VectorParameter <- function(o) { inherits(o, "VectorParameter") }
#' @export
is.MatrixParameter <- function(o) { inherits(o, "MatrixParameter") }


#' @export
is.WithCustomVecParams <- function(o) { inherits(o, "_WithCustomVecParams") }


## constant
#' @export
is.Fixed <- function(o) { inherits(o, "_Fixed") || inherits(o, "_Zeros") || inherits(o, "_Ones") || inherits(o, "_Identity") }
#' @export
is.Zeros <- function(o) { inherits(o, "_Zeros") }
#' @export
is.Ones <- function(o) { inherits(o, "_Ones") }
#' @export
is.Identity <- function(o) { inherits(o, "_Identity") }

## Properties
#' @export
is.AllEqual <- function(o) { inherits(o, "_AllEqual") || inherits(o, "_ScalarDiagonal") || inherits(o, "ScalarParameter") }
#' @export
is.NonNegative <- function(o) { inherits(o, "_NonNegative") }
#' @export
is.Diagonal <- function(o) { inherits(o, "_Diagonal") || inherits(o, "_ScalarDiagonal") }
#' @export
is.ScalarDiagonal <- function(o) { inherits(o, "_ScalarDiagonal") }
#' @export
is.Symmetric <- function(o) { inherits(o, "_Symmetric") }
# upper triangular excluding diagonal
#' @export
is.UpperTriangular <- function(o) { inherits(o, "_UpperTriangular") }
# upper triangular including diagonal
#' @export
is.UpperTriangularWithDiagonal <- function(o) { inherits(o, "_UpperTriangularWithDiagonal") || inherits(o, "_CholeskiFactor") }
#' @export
is.WithNonNegativeDiagonal <- function(o) { inherits(o, "_WithNonNegativeDiagonal") || inherits(o, "_CholeskiFactor") }
# lower triangular excluding diagonal
#' @export
is.LowerTriangular <- function(o) { inherits(o, "_LowerTriangular") }
# lower triangular with diagonal
#' @export
is.LowerTriangularWithDiagonal <- function(o) { inherits(o, "_LowerTriangularWithDiagonal") }
#' @export
is.Omitted <- function(o) { inherits(o, "_Omitted") }

#' @export
is.CholeskiFactor <- function(o) { inherits(o, "_CholeskiFactor") }

#' @export
is.Schur <- function(o) { inherits(o, "_Schur") }

#' @export
is.Transformable <- function(o) { inherits(o, "_Transformable") }

#' @export
is.Transformed <- function(o) { inherits(o, "_Transformed") }

#' @export
is.SemiPositiveDefinite <- function(o) { inherits(o, "_SemiPositiveDefinite") }

### .. Overwrite `base::diag<-` for the case of k=1
`diag<-` <- function (x, value)
{
  dx <- dim(x)
  if( is.null(dx) ) {
    # x is a vector, assign the first element
    if(length(value) != 1) {
      stop("replacement diagonal has wrong length")
    } else {
      x[1] <- value
    }
  } else if (length(dx) != 2L) {
    stop("only matrix diagonals can be replaced")
  } else {
    len.i <- min(dx)
    len.v <- length(value)
    if (len.v != 1L && len.v != len.i)
      stop("replacement diagonal has wrong length")
    if (len.i) {
      i <- seq_len(len.i)
      x[cbind(i, i)] <- value
    }
  }
  x
}

### ..Load a from b -----------
`%load%` <- `<-`

### ..Store a to b ------------
`%store%` <- function(a,b) eval(substitute(b<-a), parent.frame())

#' Load a PCM parameter from a vector of the variable parameters in a model.
#'
#' @details This S3 generic function has both, a returned value and side effects.
#' @return an integer equaling the number of elemnents read from vecParams.
#' In the case of type=="custom", the number of indices bigger than offset returned by the function indices(offset, k).
#' @export
PCMParamLoadOrStore <- function(о, vecParams, offset, k, R, load, parentModel = NULL) {
  UseMethod("PCMParamLoadOrStore", о)
}

#' @export
PCMParamLoadOrStore.ScalarParameter <- function(о, vecParams, offset, k, R, load, parentModel = NULL) {
  `%op%` <- if(load) `%load%` else `%store%`
  if(is.Fixed(o) || is.Omitted(o)) {
    # do nothing
    num <- 0L
  } else if(is.WithCustomVecParams(o)) {
    indices <- attr(o, "indices")
    if(is.function(indices)) {

      if(is.Global(o)) {
        ind <- indices(offset, k)
        eval(substitute(o[] %op% vecParams[ind]), parent.frame())
        num <- length(unique(ind[ind > offset]))
      } else {
        num <- 0L
        for(r in 1:R) {
          ind <- indices(offset + num, k)
          eval(substitute(o[r] %op% vecParams[ind]), parent.frame())
          num <- num + length(unique(ind[ind > offset + num]))
        }
      }
    } else {
      stop("ERR:02701:PCMBase:PCMParam.R:PCMParamLoadOrStore.ScalarParameter:: indices should be a function(offset, k) returning an integer.")
    }
  } else {
    if(is.Global(o)) {
      eval(substitute(o[] %op% vecParams[offset + num]), parent.frame())
      num <- 1L
    } else {
      eval(substitute(o[1:R] %op% vecParams[offset + (1:R)]), parent.frame())
      num <- R
    }
  }
  num
}

#' @export
PCMParamLoadOrStore.VectorParameter <- function(o, vecParams, offset, k, R, load, parentModel = NULL) {
  `%op%` <- if(load) `%load%` else `%store%`

  if(is.Fixed(o) || is.Omitted(o)) {
    num <- 0
  } else if(is.WithCustomVecParams(o)) {
    indices <- attr(o, "indices", exact = TRUE)
    mask <- attr(o, "mask", exact=TRUE)
    if(is.function(indices)) {

      if(is.Global(o)) {
        ind <- indices(offset, k)
        eval(substitute(o[mask] %op% vecParams[ind]), parent.frame())
        num <- length(unique(ind[ind > offset]))
      } else {
        num <- 0L
        for(r in 1:R) {
          ind <- indices(offset + num, k)
          eval(substitute(o[mask, r] %op% vecParams[ind]), parent.frame())
          num <- num + length(unique(ind[ind > offset + num]))
        }
      }
    } else {
      stop("ERR:02711:PCMBase:PCMParam.R:PCMParamLoadOrStore.VectorParameter:: indices should be a function(offset, k) returning an integer vector.")
    }
  } else if(is.AllEqual(o)) {
    if(load) {
      mask <- rep(TRUE, k)
    } else {
      mask <- rep(FALSE, k)
      mask[1L] <- TRUE
    }

    if(is.Global(o)) {
      eval(substitute(o[mask] %op% vecParams[offset + 1L]), parent.frame())
      num <- 1L
    } else {
      num <- 0L
      for(r in 1:R) {
        eval(substitute(o[mask, r] %op% vecParams[offset + num + 1L]), parent.frame())
        num <- num + 1L
      }
    }
  } else {
    # type = Full

    mask <- rep(TRUE, k)
    if(is.Global(o)) {
      eval(substitute(o[mask] %op% vecParams[offset + (1:k)]), parent.frame())
      num <- k
    } else {
      num <- 0L
      for(r in 1:R) {
        eval(substitute(o[, r][mask] %op% vecParams[offset + num + (1:k)]), parent.frame())
        num <- num + k
      }
    }
  }
  num
}

#' @export
PCMParamLoadOrStore.MatrixParameter <- function(o, vecParams, offset, k, R, load, parentModel = NULL) {
  `%op%` <- if(load) `%load%` else `%store%`

  if(is.Fixed(o) || is.Omitted(o)) {
    # nothing to do
    num <- 0L
  } else if(is.WithCustomVecParams(o)) {
    indices <- attr(o, "indices", exact = TRUE)
    mask <- attr(o, "mask", exact=TRUE)
    if(is.function(indices)) {
      if(is.Global(o)) {
        ind <- indices(offset, k)
        eval(substitute(o[mask] %op% vecParams[ind]), parent.frame())
        num <- length(unique(ind[ind > offset]))
      } else {
        num <- 0L
        for(r in 1:R) {
          ind <- indices(offset + num, k)
          eval(substitute(o[,,r][mask] %op% vecParams[ind]), parent.frame())
          num <- num + length(unique(ind[ind > offset + num]))
        }
      }
    } else {
      stop("ERR:02721:PCMBase:PCMParam.R:PCMParamLoadOrStore.MatrixParameter:: indices should be a function(offset, k) returning an integer vector.")
    }
  } else if(is.ScalarDiagonal(o)) {
    if(load) {
      mask <- diag(TRUE, k, k)
    } else {
      mask <- matrix(FALSE, k, k)
      mask[1L, 1L] <- TRUE
    }
    if(is.Global(o)) {
      eval(substitute(o[mask] %op% vecParams[offset + 1L]), parent.frame())
      num <- 1L
    } else {
      num <- 0L
      for(r in 1:R) {
        eval(substitute(o[, , r][mask] %op% vecParams[offset + num + 1L]), parent.frame())
        num <- num + 1L
      }
    }
  } else if(is.Diagonal(o)) {
    mask <- diag(TRUE, k, k)
    if(is.Global(o)) {
      eval(substitute(o[mask] %op% vecParams[offset + (1:k)]), parent.frame())
      num <- k
    } else {
      num <- 0L
      for(r in 1:R) {
        eval(substitute(o[,, r][mask] %op% vecParams[offset + num + (1:k)]), parent.frame())
        num <- num + k
      }
    }
  } else if(is.UpperTriangular(o)) {
    # without diagonal
    numForOneMatrix <- k*(k-1)/2
    if(is.Global(o)) {
      mask <- upper.tri(o)
      eval(substitute(o[mask] %op% vecParams[offset + (1:numForOneMatrix)]), parent.frame())
      num <- numForOneMatrix
    } else {
      mask <- upper.tri(o[,,1L])
      num <- 0L
      for(r in 1:R) {
        eval(substitute(o[,,r][mask] %op% vecParams[offset + num + (1:numForOneMatrix)]), parent.frame())
        num <- num + numForOneMatrix
      }
    }
  } else if(is.UpperTriangularWithDiagonal(o)) {
    numForOneMatrix <- k*(k+1)/2
    if(is.Global(o)) {
      mask <- upper.tri(o, diag = TRUE)
      eval(substitute(o[mask] %op% vecParams[offset + (1:numForOneMatrix)]), parent.frame())
      num <- numForOneMatrix
    } else {
      mask <- upper.tri(o[,,1L], diag = TRUE)
      num <- 0L
      for(r in 1:R) {
        eval(substitute(o[,,r][mask] %op% vecParams[offset + num + (1:numForOneMatrix)]), parent.frame())
        num <- num + numForOneMatrix
      }
    }
  } else if(is.LowerTriangular(o)) {
    # without diagonal
    numForOneMatrix <- k*(k-1)/2
    if(is.Global(o)) {
      mask <- lower.tri(o)
      eval(substitute(o[mask] %op% vecParams[offset + (1:numForOneMatrix)]), parent.frame())
      num <- numForOneMatrix
    } else {
      mask <- lower.tri(o[,,1L])
      num <- 0L
      for(r in 1:R) {
        eval(substitute(o[,,r][mask] %op% vecParams[offset + num + (1:numForOneMatrix)]), parent.frame())
        num <- num + numForOneMatrix
      }
    }
  } else if(is.LowerTriangularWithDiagonal(o)) {
    numForOneMatrix <- k*(k+1)/2
    if(is.Global(o)) {
      mask <- lower.tri(o, diag = TRUE)
      eval(substitute(o[mask] %op% vecParams[offset + (1:numForOneMatrix)]), parent.frame())
      num <- numForOneMatrix
    } else {
      mask <- lower.tri(o[,,1L], diag = TRUE)
      num <- 0L
      for(r in 1:R) {
        eval(substitute(o[,,r][mask] %op% vecParams[offset + num + (1:numForOneMatrix)]), parent.frame())
        num <- num + numForOneMatrix
      }
    }
  } else if(is.Symmetric(o)) {
    numForOneMatrix <- k*(k+1)/2
    if(is.Global(o)) {
      mask <- upper.tri(o, diag = TRUE)
      maskLower <- lower.tri(o)
      eval(substitute(о[mask] %op% vecParams[offset + (1:numForOneMatrix)]), parent.frame())
      if(load) {
        eval(substitute({
          o[maskLower] <- 0.0
          o <- o + t(o)
          diag(o) <- 0.5 * diag(o)
        }), parent.frame())
      }
      num <- numForOneMatrix
    } else {
      mask <- upper.tri(o[,,1L], diag = TRUE)
      maskLower <- lower.tri(o[,,1L])
      num <- 0L
      for(r in 1:R) {
        eval(substitute(o[,,r][mask] %op% vecParams[offset + num + (1:numForOneMatrix)]), parent.frame())
        if(load) {
          eval(substitute({
            o[,,r][maskLower] <- 0.0
            o[,,r] <- o[,,r] + t(o[,,r])
            diag(o[,,r]) <- 0.5 * diag(o[,,r])
          }), parent.frame())
        }
        num <- num + numForOneMatrix
      }
    }
  } else {
    # Full matrix
    numForOneMatrix <- k * k
    mask <- matrix(TRUE, k, k)

    if(is.Global(o)) {
      eval(substitute(o[mask] %op% vecParams[offset + (1:numForOneMatrix)]), parent.frame())
      num <- numForOneMatrix
    } else {
      num <- 0L
      for(r in 1:R) {
        eval(substitute(o[,,r][mask] %op% vecParams[offset +  num + (1:numForOneMatrix)]), parent.frame())
        num <- num + numForOneMatrix
      }
    }
  }
  num
}

#' @export
PCMParamLoadOrStore.PCM <- function(o, vecParams, offset, k, R, load, parentModel = NULL) {
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))

  p <- 0

  if(is.Fixed(o) || is.Omitted(o)) {
    # do nothing
  } else {
    for(name in names(o)) {
      if(is.Global(o[[name]]) && !is.null(parentModel)) {
        # this is a nested model of parentModel, hence
        # o[[name]] must have been loaded by the parentModel, because all global
        # parameters by convention should precede the other entries.
        if(load) {
          eval(substitute( o[[name]][] <- parentModel[[name]][] ), parent.frame() )
          # no need to update p
        } else {
          # nothing to do because the parentModel should already have taken care
          # to store the parameter into vecParams.
        }
      } else {
        # naming the first parameter (o) seems to fail:
        #
        # p <- p + eval(substitute(PCMParamLoadOrStore(
        #   o = o[[name]],
        #   vecParams = vecParams,
        #   offset = offset + p,
        #   k = k,
        #   R = R,
        #   load = load)), parent.frame())
        #
        p <- p + eval(substitute(PCMParamLoadOrStore(
          o[[name]],
          vecParams,
          offset + p,
          k,
          R,
          load,
          o)), parent.frame())
      }

    }
  }
  p
}



#' Count the number of free parameters associated with a PCM or a PCM-parameter
#' @param o a PCM model object or a parameter of a PCM object
#' @param countRegimeChanges logical indicating if regime changes should be counted.
#' If TRUE, the default implementation would add \code{PCMNumRegimes(model) - 1}.
#' Default FALSE.
#' @param countModelTypes logical indicating whether the model type should be
#'  counted. If TRUE the default implementation will add +1 only if there are more than
#'  one modelTypes (\code{length(attr(model, "modelTypes", exact = TRUE)) > 1}),
#'  assuming that all regimes are regimes of the same model type (e.g. OU). The
#'  implementation for MRG models will add +1 for every regime if there are more than
#'  one modelTypes. Default FALSE.
#' @return an integer
#' @export
PCMParamCount <- function(o, countRegimeChanges = FALSE, countModelTypes = FALSE, offset = 0, k = 1, R = 1, parentModel = NULL) {
  UseMethod("PCMParamCount", o)
}

#' @export
PCMParamCount.ScalarParameter <- function(o, countRegimeChanges = FALSE, countModelTypes = FALSE,  offset = 0, k = 1, R = 1, parentModel = NULL) {
  if(is.Fixed(o)) {
    0L
  } else if(is.WithCustomVecParams(o)) {
    indices <- attr(o, "indices", exact = TRUE)
    p0 <- offset
    if(is.Global(o)) {
      ind <- indices(offset, k)
      offset <- offset + length(unique(ind[ind > offset]))
    } else {
      for(r in 1:R) {
        ind <- indices(offset, k)
        offset <- offset + length(unique(ind[ind > offset]))
      }
    }
    offset - p0
  } else {
    if(is.Global(o)) {
      1L
    } else {
      R
    }
  }
}

#' @export
PCMParamCount.VectorParameter <- function(o, countRegimeChanges = FALSE, countModelTypes = FALSE, offset = 0, k = 1, R = 1, parentModel = NULL) {
  if(is.Fixed(o) || is.Omitted(o)) {
    0L
  } else if(is.AllEqual(o)) {
    if(is.Global(o)) {
      1L
    } else {
      R
    }
  } else if(is.WithCustomVecParams(o)) {
    indices <- attr(o, "indices", exact = TRUE)
    if(is.function(indices)) {
      p0 <- offset
      if(is.Global(o)) {
        ind <- indices(offset, k)
        offset <- offset + length(unique(ind[ind > offset]))
      } else {
        for(r in 1:R) {
          ind <- indices(offset, k)
          offset <- offset + length(unique(ind[ind > offset]))
        }
      }
      offset - p0
    } else {
      stop("ERR:02731:PCMBase:PCMParam.R:PCMParamCount.VectorParameter:: indices should be a function(offset, k) returning an integer vector.")
    }
  } else {
    if(is.Global(o)) {
      k
    } else {
      R * k
    }
  }
}

#' @export
PCMParamCount.MatrixParameter <- function(o, countRegimeChanges = FALSE, countModelTypes = FALSE,  offset = 0, k = 1, R = 1, parentModel = NULL) {
  if(is.Fixed(o) || is.Omitted(o)) {
    # nothing to do
    0L
  } else if(is.WithCustomVecParams(o)) {
    indices <- attr(o, "indices", exact = TRUE)
    if(is.function(indices) ) {
      p0 <- offset
      if(is.Global(o)) {
        ind <- indices(offset, k)
        offset <- offset + length(unique(ind[ind > offset]))
      } else {
        for(r in 1:R) {
          ind <- indices(offset, k)
          offset <- offset + length(unique(ind[ind > offset]))
        }
      }
      offset - p0
    } else {
      stop("ERR:02741:PCMBase:PCMParam.R:PCMParamCount.MatrixParameter:: indices should be a function(offset, k) returning an integer vector.")
    }
  } else if(is.ScalarDiagonal(o)) {
    if(is.Global(o)) {
      1L
    } else {
      R
    }
  } else if(is.Diagonal(o)) {
    if(is.Global(o)) {
      k
    } else {
      R * k
    }
  } else if(is.UpperTriangular(o)) {
    # without diagonal
    if(is.Global(o)) {
      k*(k-1)/2
    } else {
      R * k*(k-1)/2
    }
  } else if(is.UpperTriangularWithDiagonal(o)) {
    # with diagonal
    if(is.Global(o)) {
      k*(k+1)/2
    } else {
      R * k*(k+1)/2
    }
  } else if(is.LowerTriangular(o)) {
    # without diagonal
    if(is.Global(o)) {
      k*(k-1)/2
    } else {
      R * k*(k-1)/2
    }
  } else if(is.LowerTriangularWithDiagonal(o)) {
    if(is.Global(o)) {
      k*(k+1)/2
    } else {
      R * k*(k+1)/2
    }
  } else if(is.Symmetric(o)) {
    if(is.Global(o)) {
      k*(k+1)/2
    } else {
      R * k*(k+1)/2
    }
  } else {
    # Full matrix
    if(is.Global(o)) {
      k*k
    } else {
      R * k*k
    }
  }
}

#' @export
PCMParamCount.PCM <- function(o, countRegimeChanges = FALSE, countModelTypes = FALSE, offset = 0, k = 1, R = 1, parentModel = NULL) {
  k <- attr(o, "k")
  R <- length(attr(o, "regimes"))

  p0 <- offset
  if(is.Fixed(o) || is.Omitted(o)) {
    # do nothing
  } else {
    for(name in names(o)) {
      if(is.Global(o[[name]]) && !is.null(parentModel)) {
        # Global parameters are already counted in the parentModel, so don't recount them
      } else {
        offset <- offset + PCMParamCount(o[[name]], offset = offset, k = k, R = R, parentModel = o)
      }
    }
    if(countRegimeChanges) {
      # we don't count the root as a parameters, tha'ts why we substract one.
      offset <- offset + PCMNumRegimes.PCM(o) - 1L
    }
    if(countModelTypes) {
      # assume that all regimes have the same model-type. If there is only one
      # model type than this is not counted as a parameter.
      if(length(attr(o, "modelTypes", exact = TRUE)) > 1L) {
        offset <- offset + 1L
      }
    }
  }
  unname(offset - p0)
}


#' @export
PCMParamGetShortVector <- function(o, k = 1, R = 1, ...) {
  UseMethod("PCMParamGetShortVector", o)
}

#' @export
PCMParamGetShortVector.PCM <- function(o, k, R, ...) {
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))

  vec <- double(PCMParamCount(o, FALSE, FALSE, 0L, k, R))
  PCMParamLoadOrStore(o, vecParams = vec, offset = 0L, k = k, R = R, load = FALSE)
  vec
}

#' @export
PCMParamGetShortVector.default <- function(o, k, R, ...) {
  vec <- double(PCMParamCount(o, FALSE, FALSE, 0L, k, R))
  PCMParamLoadOrStore(o, vecParams = vec, offset = 0L, k = k, R = R, load = FALSE)
  vec
}

#' Set model parameters from a named list
#' @param tree a phylo object (possible future use)
#' @param model a PCM model object
#' @param params a named list with elements among the names found in model
#' @param inplace logical indicating if the parameters should be set "inplace" for the
#' model object in the calling environment or a new model object with the parameters set
#' as specified should be returned. Defaults to TRUE
#' @return If inplace is TRUE, the function only has a side effect of setting the
#' parameters of the model object in the calling environment; otherwise the function
#' returns a modified copy of the model object.
#' @export
PCMParamSetByName <- function(model, params, inplace = TRUE, replaceWholeParameters = FALSE, ...) {
  UseMethod("PCMParamSetByName", model)
}

#' @export
PCMParamSetByName.PCM <- function(model, params, inplace = TRUE, replaceWholeParameters = FALSE, ...) {
  for(name in names(params)) {
    if(! (name %in% names(model)) ) {
      stop(paste0("ERR:02751:PCMBase:PCMParam.R:PCMParamSetByName.PCM:: ", name,
                  " is not a settable parameter of the model."))
    }

    if(is.PCM(model[[name]])) {
      if(! is.PCM(params[[name]]) ) {
        stop(paste0("ERR:02752:PCMBase:PCMParam.R:PCMParamSetByName.PCM::model[['", name, "']] is a nested PCM object but params[['", name, "']] is not."))
      } else if( !identical(attr(model, "k", exact = TRUE), attr(params[[name]], "k", exact = TRUE)) ) {
        stop(paste0("ERR:02753:PCMBase:PCMParam.R:PCMParamSetByName.PCM:: model[['", name, "']] has a different number of traits (k) from params[['", name, ",]]."))
      } else {
        if(inplace) {
          eval(substitute(model[[name]] <- params[[name]]), parent.frame())
        } else {
          model[[name]] <- params[[name]]
        }
      }
    } else {
      if( !identical(length(model[[name]]), length(params[[name]])) && !replaceWholeParameters ) {
        stop(paste0("ERR:02754:PCMBase:PCMParam.R:PCMParamSetByName.PCM:: params[[", name,
                    "]] is not the same length as model[[", name, "]]; ",
                    "length(params[[", name, "]])=", length(params[[name]]),
                    ", length(model[[", name, "]])=", length(model[[name]]), "."))
      }
      if( !identical(dim(model[[name]]), dim(params[[name]])) && !replaceWholeParameters ) {
        stop(paste0("ERR:02755:PCMBase:PCMParam.R:PCMParamSetByName.PCM:: params[[", name,
                    "]] is not the same dimension as model[[", name, "]]; ",
                    "dim(params[[", name, "]])=", str(dim(params[[name]])),
                    ", length(model[[", name, "]])=", str(dim(model[[name]])), "."))
      }

      if(replaceWholeParameters) {
        # This will keep the current attributes of model[[name]]
        if(inplace) {
          eval(substitute(model[[name]] <- params[[name]]), parent.frame())
        } else {
          model[[name]] <- params[[name]]
        }
      } else {
        # This will keep the current attributes of model[[name]]
        if(inplace) {
          eval(substitute(model[[name]][] <- params[[name]]), parent.frame())
        } else {
          model[[name]][] <- params[[name]]
        }
      }
    }
  }

  if(!inplace) {
    model
  }
}

#' The object representing a lower limit for a given type
#'
#' @description \code{PCMParamLowerLimit} and \code{PCMParamUpperLimit} are S3
#' generic functions.
#'
#' @param o an object such as a VectorParameter a MatrixParameter or a PCM.
#' @param k integer denoting the number of traits
#' @param R integer denoting the number of regimes in the model in which o
#' belongs to.
#' @param ... additional arguments (optional or future use).
#' @return an object of the same S3 class as o representing a lower limit for
#' the class.
#' @export
PCMParamLowerLimit <- function(o, k, R, ...) {
  UseMethod("PCMParamLowerLimit", o)
}

#' @export
PCMParamLowerLimit.PCM <- function(o, k, R,  ...) {
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))

  vecParamsLowerLimit <- attr(o, "vecParamsLowerLimit", exact = TRUE)

  if( !is.null(vecParamsLowerLimit) ) {
    PCMParamLoadOrStore(o, vecParams = vecParamsLowerLimit, offset = 0, k = k, R = R, load = TRUE)
  } else {
    PCMParamSetByName(o, params = lapply(o, PCMParamLowerLimit, k = k, R = R, ...), inplace = TRUE)
  }
  o
}

#' @export
PCMParamLowerLimit.default <- function(o, k, R,  ...) {
  vecParamsLowerLimit <- attr(o, "vecParamsLowerLimit", exact = TRUE)
  if( !is.null(vecParamsLowerLimit) ) {
    PCMParamLoadOrStore(o, vecParams = vecParamsLowerLimit, offset = 0, k = k, R = R, load = TRUE)
  } else {
    nParams <- PCMParamCount(o, countRegimeChanges = FALSE, countModelTypes = FALSE,
                             offset = 0, k = k, R = R)
    vecParamsLowerLimit <- rep(getOption("PCMBase.ParamValue.LowerLimit", -10.0), nParams)
    PCMParamLoadOrStore(o, vecParams = vecParamsLowerLimit, offset = 0, k = k, R = R, load = TRUE)
  }
  o
}

#' @export
PCMParamLowerLimit._WithNonNegativeDiagonal <- function(o, k, R, ...) {
  vecParamsLowerLimit <- attr(o, "vecParamsLowerLimit", exact = TRUE)
  if( !is.null(vecParamsLowerLimit) ) {
    PCMParamLoadOrStore(o, vecParams = vecParamsLowerLimit, offset = 0, k = k, R = R, load = TRUE)
  } else {
    o2 <- PCMParamLowerLimit.default(o, k, R, ...)

    valueLowerLimitDiagonal <- getOption("PCMBase.ParamValue.LowerLimit.NonNegativeDiagonal", 0.0)

    if(is.Global(o2)) {
      diag(o2) <- valueLowerLimitDiagonal
    } else {
      for(r in 1:R) {
        diag(o2[,,r]) <- valueLowerLimitDiagonal
      }
    }

    # at this point there is a risk that we have overwritten some fixed values on the diagonal(s) in o2.
    # We correct for this by taking the vecParams for o2 and loading it into o
    o2VecParams <- PCMParamGetShortVector(o2, k, R)
    PCMParamLoadOrStore(o, o2VecParams, 0L, k, R, TRUE)
  }

  o
}

#' @export
PCMParamUpperLimit <- function(o, k, R,  ...) {
  UseMethod("PCMParamUpperLimit", o)
}

#' @export
PCMParamUpperLimit.PCM <- function(o, k, R,  ...) {
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))

  vecParamsUpperLimit <- attr(o, "vecParamsUpperLimit", exact = TRUE)

  if( !is.null(vecParamsUpperLimit) ) {
    PCMParamLoadOrStore(o, vecParams = vecParamsUpperLimit, offset = 0, k = k, R = R, load = TRUE)
  } else {
    PCMParamSetByName(o, params = lapply(o, PCMParamUpperLimit, k = k, R = R, ...), inplace = TRUE)
  }
  o
}

#' @export
PCMParamUpperLimit.default <- function(o, k, R,  ...) {
  vecParamsUpperLimit <- attr(o, "vecParamsUpperLimit", exact = TRUE)
  if( !is.null(vecParamsUpperLimit) ) {
    PCMParamLoadOrStore(o, vecParams = vecParamsUpperLimit, offset = 0, k = k, R = R, load = TRUE)
  } else {
    nParams <- PCMParamCount(o, countRegimeChanges = FALSE, countModelTypes = FALSE,
                             offset = 0, k = k, R = R)
    vecParamsUpperLimit <- rep(getOption("PCMBase.ParamValue.UpperLimit", 10.0), nParams)
    PCMParamLoadOrStore(o, vecParams = vecParamsUpperLimit, offset = 0, k = k, R = R, load = TRUE)
  }
  o
}

#' Generate a random parameter vector for a model using uniform distribution between its lower and upper bounds.
#' @param o a PCM model object or a parameter
#' @param n an integer specifying the number of random vectors to generate
#' @param argsPCMParamLowerLimit,argsPCMParamUpperLimit named lists of arguments passed to
#' \code{\link{PCMParamLowerLimit}} and \code{\link{PCMParamUpperLimit}}.
#' @return if n = 1, a numeric vector of length \code{PCMParamCount(o)}; if n > 1,
#' a numeric matrix of dimension n x \code{PCMParamCount(o)}.
#' @seealso PCMParamUpperLimit PCMParamLowerLimit PCMParamGetShortVector
#' @export
PCMParamRandomVecParams <- function(o, k, R, n = 1L,
                                    argsPCMParamLowerLimit = NULL,
                                    argsPCMParamUpperLimit = NULL) {
  UseMethod("PCMParamRandomVecParams", o)
}

#' @importFrom stats runif
#' @export
PCMParamRandomVecParams.default <- function(o, k, R, n = 1L,
                                        argsPCMParamLowerLimit = NULL,
                                        argsPCMParamUpperLimit = NULL) {

  if(is.PCM(o)) {
    k <- attr(o, "k", exact = TRUE)
    R <- length(attr(o, "regimes", exact = TRUE))
  }

  lowerModel <- do.call(PCMParamLowerLimit, c(list(o, k = k, R = R), argsPCMParamLowerLimit))
  lowerVecParams <- PCMParamGetShortVector(lowerModel, k = k, R = R)

  upperModel <- do.call(PCMParamUpperLimit, c(list(o, k = k, R = R), argsPCMParamUpperLimit))
  upperVecParams <- PCMParamGetShortVector(upperModel, k = k, R = R)

  p <- PCMParamCount(o, k = k, R = R)
  res <- sapply(1:p, function(i) {
    runif(n, lowerVecParams[i], upperVecParams[i])
  })
  res
}

#' Generate a default object of a given PCM model type or parameter type
#' @param spec any object having a class attribute. The value of this object is not
#' used, but its class is used for method-dispatch.
#' @param model a PCM object used to extract attributes needed for creating a
#' default object of class specified in \code{class(spec)}, such as the number of
#' traits (k) or the regimes and the number of regimes;
#' @return a parameter or a PCM object.
#' @export
PCMDefaultObject <- function(spec, model, ...) {
  UseMethod("PCMDefaultObject", spec)
}

#' @export
PCMDefaultObject.PCM <- function(spec, model, ...) {
  o <- spec
  PCMParamSetByName(o, lapply(spec, PCMDefaultObject, model = o), replaceWholeParameters = TRUE, ...)
  o[] <- o[!sapply(o, is.null)]
  attr(o, "spec") <- spec
  o
}

#' @export
PCMDefaultObject.ScalarParameter <- function(spec, model, ...) {
  regimes <- attr(model, "regimes", exact = TRUE)
  R <- length(regimes)

  if(is.Omitted(spec)) {
    o <- NULL
  } else if(is.Global(spec)) {
    o <- structure(
      0.0,
      class = c(class(spec), "numeric"),
      description = attr(spec, "description", exact = TRUE),
      vecParamsLowerLimit = attr(spec, "vecParamsLowerLimit", exact = TRUE),
      vecParamsUpperLimit = attr(spec, "vecParamsUpperLimit", exact = TRUE))
  } else {
    o <- structure(
      double(R),
      names = regimes,
      class = c(class(spec), "numeric"),
      description = attr(spec, "description", exact = TRUE),
      vecParamsLowerLimit = attr(spec, "vecParamsLowerLimit", exact = TRUE),
      vecParamsUpperLimit = attr(spec, "vecParamsUpperLimit", exact = TRUE))
  }
  o
}

#' @export
PCMDefaultObject.VectorParameter <- function(spec, model, ...) {
  k <- attr(model, "k", exact = TRUE)
  regimes <- attr(model, "regimes", exact = TRUE)
  R <- length(regimes)

  if(is.Omitted(spec)) {
    o <- NULL
  } else if(is.Global(spec)) {
    o <- structure(
      double(k),
      class = c(class(spec), "numeric"),
      description = attr(spec, "description", exact = TRUE),
      vecParamsLowerLimit = attr(spec, "vecParamsLowerLimit", exact = TRUE),
      vecParamsUpperLimit = attr(spec, "vecParamsUpperLimit", exact = TRUE))
  } else {
    o <- structure(
      array(0.0, dim = c(k, R), dimnames = list(NULL, regimes)),
      class = c(class(spec), "matrix"),
      description = attr(spec, "description", exact = TRUE),
      vecParamsLowerLimit = attr(spec, "vecParamsLowerLimit", exact = TRUE),
      vecParamsUpperLimit = attr(spec, "vecParamsUpperLimit", exact = TRUE))
  }
  if(is.Ones(spec)) {
    o[] <- 1.0
  }
  o
}

#' @export
PCMDefaultObject.MatrixParameter <- function(spec, model, ...) {
  k <- attr(model, "k", exact = TRUE)
  regimes <- attr(model, "regimes", exact = TRUE)
  R <- length(regimes)

  if(is.Omitted(spec)) {
    o <- NULL
  } else if(is.Global(spec)) {
    o <- structure(
      matrix(0.0, k, k),
      class = c(class(spec), "matrix"),
      description = attr(spec, "description", exact = TRUE),
      vecParamsLowerLimit = attr(spec, "vecParamsLowerLimit", exact = TRUE),
      vecParamsUpperLimit = attr(spec, "vecParamsUpperLimit", exact = TRUE))
  } else {
    o <- structure(
      array(0.0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
      class = class(spec),
      description = attr(spec, "description", exact = TRUE),
      vecParamsLowerLimit = attr(spec, "vecParamsLowerLimit", exact = TRUE),
      vecParamsUpperLimit = attr(spec, "vecParamsUpperLimit", exact = TRUE))
  }
  if(is.Identity(spec)) {
    if(is.Global(spec)) {
      diag(o) <- 1.0
    } else {
      for(r in 1:R) {
        diag(o[,,r]) <- 1.0
      }
    }
  } else if(is.Ones(spec)) {
    o[] <- 1.0
  }
  o
}

#' @export
PCMApplyTransformation._CholeskiFactor <- function(o, ...) {
  # when assigning to o, we use o[] <-  instead of just o <- , in order to
  # preserve the attributes
  if(is.Global(o)) {
    o[] <- o %*% t(o)
  } else {
    for(r in 1:dim(o)[3]) {
      o[,,r] <- o[,,r] %*% t(o[,,r])
    }
  }
  classes <- class(o)
  classes <- classes[!(classes %in% c("_UpperTriangular", "_UpperTriangularWithDiagonal", "_CholeskiFactor", "_Transformable", "matrix"))]
  if(is.WithNonNegativeDiagonal(o)) {
    classes <- c(classes, "_SemiPositiveDefinite")
  }

  classes <- c(classes, "_Transformed")

  if(is.Global(o)) {
    # o is a matrix, so we need to set its "matrix" class explicitly
    classes <- c(classes, "matrix")
  }
  class(o) <- classes
  o
}

#' @export
PCMApplyTransformation._Schur <- function(o, ...) {
  # Assuming o is a MatrixParameter, we transform each matrix M in o as follows:
  # use the upper triangle without the diagonal of M as rotation angles for a Givens
  # rotation to obtain an orthoganal matrix Q;
  # Then, use the lower triangle with the diagonal of M as a matrix T
  # Return the product Q %*% t(T) %*% t(Q).
  # The returned matrix will have all of its eigenvalues equal to the elements
  # on the diaogonal of M (and T). If M is upper triangular, then T will be diagonal
  # and the returned matrix will be symmetric. If the elements on the diagonal
  # of M are positive, so will be the eigenvalues of the returned matrix.

  transformMatrix <- function(M) {
    k <- nrow(M)
    if(k == 1) {
      Q <- matrix(1.0, 1, 1)
    } else {
      angles <- M[upper.tri(M)]
      Q <- .par.transform.orth.matrix.givens(angles, k)$Q
    }
    M[upper.tri(M)] <- 0.0
    Q %*% t(M) %*% t(Q)
  }

  if(is.Global(o)) {
    o[] <- transformMatrix(o)
  } else {
    for(r in 1:dim(o)[3]) {
      o[,,r] <- transformMatrix(o[,,r])
    }
  }
  classes <- class(o)
  classes <- classes[!(classes %in% c("_UpperTriangular", "_UpperTriangularWithDiagonal", "_Schur", "_Transformable", "matrix"))]
  if(is.UpperTriangularWithDiagonal(o)) {
    classes <- c(classes, "_Symmetric")
  }
  if(is.WithNonNegativeDiagonal(o)) {
    classes <- c(classes, "_SemiPositiveDefinite")
  }

  classes <- c(classes, "_Transformed")

  if(is.Global(o)) {
    # o is a matrix, so we need to set its "matrix" class explicitly
    classes <- c(classes, "matrix")
  }
  class(o) <- classes
  o
}
