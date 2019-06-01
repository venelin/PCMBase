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
#'
#' @return an object of S3 class PCMTable
#' @details This is an S3 generic.
#'
#' @export
PCMTable <- function(model, skipGlobalRegime = FALSE) {
  if(is.PCM(model)) {
    UseMethod("PCMTable", model)
  } else {
    stop("PCMTable:: model should be a PCM object.")
  }
}

TableForRegime <- function(model, r) {
  table <- data.table(regime = if(is.na(r)) ':global:' else as.character(r))
  for(name in names(model)) {
    if(!is.PCM(model[[name]]) && is.Global(model[[name]]) && r == ':global:') {
      value <- model[[name]]
      if(name %in% c("Sigma_x", "Sigmae_x")) {
        value <- value %*% t(value)
        name <- gsub('_x', '', name, fixed = TRUE)
      }
      table[[name]] <- list(value)
    } else if(!is.PCM(model[[name]]) && is.Local(model[[name]]) && r != ':global:') {
      value <- model[[name]]
      if(is.ScalarParameter(model[[name]])) {
        value <- model[[name]][as.character(r)]
      } else if(is.VectorParameter(model[[name]])) {
        value <- model[[name]][,as.character(r)]
      } else {
        value <- model[[name]][,,as.character(r)]
      }
      if(name %in% c("Sigma_x", "Sigmae_x")) {
        value <- value %*% t(value)
        name <- gsub('_x', '', name, fixed = TRUE)
      }
      table[[name]] <- list(value)
    }
  }
  table
}

#' @importFrom data.table rbindlist data.table setattr
#' @export
PCMTable.PCM <- function(model, skipGlobalRegime = FALSE) {
  regimes <- if(skipGlobalRegime) {
    as.character(PCMRegimes(model))
  } else {
    c(':global:', as.character(PCMRegimes(model)))
  }
  res <- rbindlist(
    lapply(regimes, function(r) {
      TableForRegime(model, r)
    }), fill=TRUE, use.names = TRUE)

  setattr(res, "model", model)
  setattr(res, "class", c("PCMTable", class(res)))
  res[, .SD, keyby=list(regime)]
  res
}

#' @importFrom data.table setcolorder
#' @export
PCMTable.MixedGaussian <- function(model) {
  res <- rbindlist(
    lapply(c(':global:', as.character(PCMRegimes(model))), function(r) {
      if(r == ':global:') {
        TableForRegime(model, r)
      } else {
        table <- PCMTable(model[[as.character(r)]], skipGlobalRegime = TRUE)
        table[, regime:=as.character(r)]
        subModelType <- match(
          class(model[[as.character(r)]])[1], attr(model, "modelTypes"))
        table[, typeId:=subModelType]
      }
    }), fill = TRUE, use.names = TRUE)

  if(!is.null(names(attr(model, "modelTypes")))) {
    res[, type:=names(attr(model, "modelTypes"))[typeId]]
    res[, typeId:=NULL]
  } else if(getOption("PCMBase.UseLettersAsTypeNames", TRUE)) {
    res[, type:=LETTERS[typeId]]
    res[, typeId:=NULL]
  } else {
    res[, type:=typeId]
    res[, typeId:=NULL]
  }
  setattr(res, "model", model)
  setattr(res, "class", c("PCMTable", class(res)))
  res[, .SD, keyby = list(regime, type, )]
  res
}

#' Latex represenataion of an object.
#' @param x an R object. Currently, objects of S3 classes MatrixParameter,
#' VectorParameter, ScalarParameter and PCMTable are supported.
#' @param ... additional arguments passed from calling functions.
#
#' @export
FormatAsLatex <- function(x, ...) {
  UseMethod("FormatAsLatex", x)
}

#' @importFrom xtable xtable
#' @export
FormatAsLatex.default <- function(x, ...) {
  if(is.null(x)) {
    ' '
  } else if(is.character(x)) {
    paste0(' ', x,' ')
  } else {
    mat <- xtable(
      as.matrix(x), align=rep("", ncol(as.matrix(x)) + 1),
      digits = getOption("digits", default = 2))
    latex <- capture.output(
      print(mat, floating=FALSE, tabular.environment="bmatrix",
            hline.after=NULL, include.rownames=FALSE,
            include.colnames=FALSE, comment = FALSE))
    paste0(' $', do.call(paste0, as.list(latex)), '$ ')
  }
}

#' @method FormatAsLatex PCMTable
#' @importFrom data.table set
#' @export
FormatAsLatex.PCMTable <- function(x, ...) {
  for(name in names(x)) {
    set(x, j = name, value = lapply(x[[name]], FormatAsLatex))
  }

  do.call(
    paste,
    as.list(
      capture.output(
        print(xtable(x),
              include.rownames= FALSE,
              sanitize.text.function=identity,
              comment = FALSE))))
}




#' @method print PCMTable
#' @export
print.PCMTable <- function(x, ...) {
  argList <- list(...)
  if(!is.null(argList$latex) && argList$latex == TRUE) {
    cat(FormatAsLatex(x, ...), "\n")
  } else {
    NextMethod()
  }
}




