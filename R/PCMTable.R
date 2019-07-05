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

#' A data.table representation of a PCM object
#' @param model a PCM object.
#' @param skipGlobalRegime logical indicating whether a raw in the returned
#' table for the global-scope parameters should be omitted (this is mostly for
#' internal use). Default (FALSE).
#' @param addTransformed logical. If TRUE (the default), columns for the
#' transformed version of the transformable parameters will be added.
#' @param removeUntransformed logical If TRUE (default), columns for the
#' untransformed version of the transformable parameters will be omitted.
#' @return an object of S3 class PCMTable
#' @details This is an S3 generic.
#'
#' @export
PCMTable <- function(
  model, skipGlobalRegime = FALSE,
  addTransformed = TRUE,
  removeUntransformed = TRUE) {

  if(is.PCM(model)) {
    UseMethod("PCMTable", model)
  } else {
    stop("PCMTable:: model should be a PCM object.")
  }
}

TableForRegime <- function(
  model, r, addTransformed = TRUE, removeUntransformed = TRUE) {

  table <- data.table(`regime` = if(is.na(r)) ':global:' else as.character(r))
  for(name in names(model)) {
    if(!is.PCM(model[[name]]) && is.Global(model[[name]]) && r == ':global:') {
      value <- model[[name]]
      nameSuffix <- ""
      if(is.Transformable(value)) {
        if(is.Schur(value)) {
          nameSuffix <- "_S"
        } else if(is.CholeskyFactor(value)) {
          nameSuffix <- "_u"
        }
      }
      table[[paste0(name, nameSuffix)]] <- list(value)
      if(addTransformed) {
        if(name %in% c("Sigma_x", "Sigmae_x", "Sigmaj_x")) {
          value <- value %*% t(value)
          name <- gsub('_x', '', name, fixed = TRUE)
          table[[name]] <- list(value)
        } else if(is.Transformable(value)) {
          value <- PCMApplyTransformation(value)
          table[[name]] <- list(value)
        }
      }
    } else if(
      !is.PCM(model[[name]]) && is.Local(model[[name]]) && r != ':global:') {

      if(is.ScalarParameter(model[[name]])) {
        value <- model[[name]][as.character(r)]
      } else if(is.VectorParameter(model[[name]])) {
        value <- model[[name]][,as.character(r)]
      } else {
        value <- model[[name]][,,as.character(r)]
      }
      nameSuffix <- ""
      if(is.Transformable(model[[name]])) {
        if(is.Schur(model[[name]])) {
          nameSuffix <- "_S"
        } else if(is.CholeskyFactor(model[[name]])) {
          nameSuffix <- "_u"
        }
      }
      table[[paste0(name, nameSuffix)]] <- list(value)

      if(addTransformed) {
        if(name %in% c("Sigma_x", "Sigmae_x", "Sigmaj_x")) {
          if(getOption("PCMBase.Transpose.Sigma_x", FALSE)) {
            value <- t(value)
          }
          value <- value %*% t(value)
          name <- gsub('_x', '', name, fixed = TRUE)
          table[[name]] <- list(value)
        } else if(is.Transformable(model[[name]])) {
          if(is.ScalarParameter(model[[name]])) {
            value <- PCMApplyTransformation(model[[name]])[as.character(r)]
          } else if(is.VectorParameter(model[[name]])) {
            value <- PCMApplyTransformation(model[[name]])[,as.character(r)]
          } else {
            value <- PCMApplyTransformation(model[[name]])[,,as.character(r)]
          }
          table[[name]] <- list(value)
        }
      }
    }
  }

  if(removeUntransformed) {
    toRemove <- sapply(names(table), function(n) {
      endsWith(n, "_x") || endsWith(n, "x}") ||
        endsWith(n, "_u")  || endsWith(n, "u}") ||
        endsWith(n, "S")
    })
    namesToRetain <- setdiff(names(table), names(table)[toRemove])
    table <- table[, namesToRetain, with=FALSE]
  }
  table
}

#' @importFrom data.table rbindlist data.table setattr
#' @export
PCMTable.PCM <- function(
  model, skipGlobalRegime = FALSE,
  addTransformed = TRUE,
  removeUntransformed = TRUE) {
  regimes <- if(skipGlobalRegime) {
    as.character(PCMRegimes(model))
  } else {
    c(':global:', as.character(PCMRegimes(model)))
  }
  res <- rbindlist(
    lapply(regimes, function(r) {
      TableForRegime(model, r, addTransformed, removeUntransformed)
    }), fill=TRUE, use.names = TRUE)

  setattr(res, "model", model)
  setattr(res, "class", c("PCMTable", class(res)))
  res
}

#' @importFrom data.table setcolorder
#' @export
PCMTable.MixedGaussian <- function(
  model, skipGlobalRegime = FALSE,
  addTransformed = TRUE,
  removeUntransformed = TRUE) {

  # NEEDED to avoid no visible binding notes during check.
  regime <- type <- typeId <- NULL
  res <- rbindlist(
    lapply(c(':global:', as.character(PCMRegimes(model))), function(r) {
      if(r == ':global:') {
        TableForRegime(model, r, addTransformed, removeUntransformed)
      } else {
        table <- PCMTable(model[[as.character(r)]], skipGlobalRegime = TRUE, addTransformed, removeUntransformed)
        table[, `regime`:=as.character(r)]
        subModelType <- match(
          class(model[[as.character(r)]])[1], attr(model, "modelTypes"))
        table[,typeId:=subModelType]
      }
    }), fill = TRUE, use.names = TRUE)

  if(!is.null(names(attr(model, "modelTypes")))) {
    res[, `type`:=names(attr(model, "modelTypes"))[typeId]]
    res[,typeId:=NULL]
  } else if(getOption("PCMBase.UseLettersAsTypeNames", TRUE)) {
    res[, `type`:=LETTERS[typeId]]
    res[,typeId:=NULL]
  } else {
    res[, `type`:=typeId]
    res[,typeId:=NULL]
  }
  setattr(res, "model", model)
  setattr(res, "class", c("PCMTable", class(res)))
  setcolorder(res, neworder = c(1, ncol(res), seq(2, ncol(res) - 1)))
  res
}

#' @method print PCMTable
#' @export
print.PCMTable <- function(x, ...) {
  argList <- list(...)
  if(!is.null(argList$xtable) && argList$xtable == TRUE) {
    cat(do.call(
      FormatTableAsLatex, c(list(x), argList[-match('xtable', names(argList))])),
      "\n")
  } else {
    NextMethod()
  }
}




