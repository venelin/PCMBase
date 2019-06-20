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
#' @return an object of S3 class PCMTable
#' @details This is an S3 generic.
#'
#' @export
PCMTable <- function(model, skipGlobalRegime = FALSE, addTransformed = TRUE) {
  if(is.PCM(model)) {
    UseMethod("PCMTable", model)
  } else {
    stop("PCMTable:: model should be a PCM object.")
  }
}

TableForRegime <- function(model, r, addTransformed = TRUE) {
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
          nameSuffix <- "_C"
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
  table
}

#' @importFrom data.table rbindlist data.table setattr
#' @export
PCMTable.PCM <- function(model, skipGlobalRegime = FALSE, addTransformed = TRUE) {
  regimes <- if(skipGlobalRegime) {
    as.character(PCMRegimes(model))
  } else {
    c(':global:', as.character(PCMRegimes(model)))
  }
  res <- rbindlist(
    lapply(regimes, function(r) {
      TableForRegime(model, r, addTransformed)
    }), fill=TRUE, use.names = TRUE)

  setattr(res, "model", model)
  setattr(res, "class", c("PCMTable", class(res)))
  res
}

#' @importFrom data.table setcolorder
#' @export
PCMTable.MixedGaussian <- function(model, skipGlobalRegime = FALSE, addTransformed = TRUE) {
  # NEEDED to avoid no visible binding notes during check.
  regime <- type <- typeId <- NULL
  res <- rbindlist(
    lapply(c(':global:', as.character(PCMRegimes(model))), function(r) {
      if(r == ':global:') {
        TableForRegime(model, r, addTransformed)
      } else {
        table <- PCMTable(model[[as.character(r)]], skipGlobalRegime = TRUE, addTransformed)
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

#' Latex represenataion of a model parameter or other found in a PCMTable object
#' @param x an R object. Currently, objects of S3 classes MatrixParameter,
#' VectorParameter, ScalarParameter and PCMTable are supported.
#' @return a character string
FormatCellAsLatex <- function(x) {
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

#' @importFrom data.table set
#' @importFrom xtable xtable
FormatPCMTableAsLatex <- function(x, argsXtable = list(), ...) {
  reGreek <- "(alpha|Alpha|beta|Beta|gamma|Gamma|delta|Delta|epsilon|Epsilon|zeta|Zeta|eta|Eta|theta|Theta|iota|Iota|kappa|Kappa|lambda|Lambda|mu|Mu|nu|Nu|xi|Xi|omicron|Omicron|pi|Pi|rho|Rho|sigma|Sigma|tau|Tau|upsilon|Upsilon|phi|Phi|chi|Chi|psi|Psi|omega|Omega)"

  for(name in names(x)) {
    set(x, j = name, value = lapply(x[[name]], FormatCellAsLatex))
  }
  setnames(x, old = names(x), new = sapply(names(x), function(n) {
    nLatexGreek <- gsub(
      '___BACKSLASH___', "\\",
      gsub(reGreek, "___BACKSLASH___\\1", n), fixed = TRUE)
    nLatexGreek <- gsub('Sigma_x', 'Sigma_{u}', nLatexGreek, fixed = TRUE)
    nLatexGreek <- gsub('Sigmae_x', 'Sigma_{e,u}', nLatexGreek, fixed = TRUE)
    nLatexGreek <- gsub('Sigmae', 'Sigma_{e}', nLatexGreek, fixed = TRUE)
    nLatexGreek <- gsub('Sigmaj_x', 'Sigma_{j,u}', nLatexGreek, fixed = TRUE)
    if(! (n %in% c('regime', 'type')) ) {
      nLatexGreek <- paste0('$',nLatexGreek,'$')
    }
    nLatexGreek
  }))

  addtorow <- list()
  addtorow$pos <- as.list(seq(1, nrow(x) - 1L))
  addtorow$command <- rep('\\vspace{2mm}  \n', length(addtorow$pos))
  do.call(
    paste,
    c(as.list(
      capture.output(
        print(do.call(xtable, c(list(x), argsXtable)),
              include.rownames= FALSE,
              add.to.row = addtorow,
              sanitize.text.function=identity,
              comment = FALSE,
              ...))),
      sep = '\n'))
}

#' @method print PCMTable
#' @export
print.PCMTable <- function(x, ...) {
  argList <- list(...)
  if(!is.null(argList$xtable) && argList$xtable == TRUE) {
    cat(do.call(
      FormatPCMTableAsLatex, c(list(x), argList[-match('xtable', names(argList))])),
      "\n")
  } else {
    NextMethod()
  }
}




