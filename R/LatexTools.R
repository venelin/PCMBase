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

#' Latex representation of a model parameter or other found in a data.table object
#' @param x an R object. Currently, character vectors of length 1,
#' numeric vectors and matrices are supported.
#'
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
            include.colnames=FALSE, comment = FALSE,
            NA.string = "."))
    paste0(' $', do.call(paste0, as.list(latex)), '$ ')
  }
}

#' Latex representation of a data.table with matrix and vectors in its cells
#' @param x a data.table
#' @param argsXtable a list (empty list by default) passed to xtable...
#' @param ... additional arguments passed to print.xtable.
#' @return a character string representing a parseable latex text.
#' @examples
#' dt <- data.table::data.table(
#'    A = list(
#'          matrix(c(2, 0, 1.2, 3), 2, 2),
#'          matrix(c(2.1, 0, 1.2, 3.2, 1.3, 3.4), 3, 2)),
#'    b = c(2.2, 3.1))
#' print(FormatTableAsLatex(dt))
#'
#' @importFrom data.table set
#' @importFrom xtable xtable
#' @export
FormatTableAsLatex <- function(x, argsXtable = list(), ...) {
  reGreek <- "(alpha|Alpha|beta|Beta|gamma|Gamma|delta|Delta|epsilon|Epsilon|zeta|Zeta|eta|Eta|theta|Theta|iota|Iota|kappa|Kappa|lambda|Lambda|mu|Mu|nu|Nu|xi|Xi|omicron|Omicron|pi|Pi|rho|Rho|sigma|Sigma|tau|Tau|upsilon|Upsilon|phi|Phi|chi|Chi|psi|Psi|omega|Omega)"

  for(name in names(x)) {
    set(x, j = name, value = lapply(x[[name]], FormatCellAsLatex))
  }
  setnames(x, old = names(x), new = sapply(names(x), function(n) {
    nLatexGreek <- gsub(
      '___BACKSLASH___', "\\",
      gsub(reGreek, "___BACKSLASH___\\1", n), fixed = TRUE)
    nLatexGreek <- gsub('_1', '^{-1}', nLatexGreek, fixed = TRUE)

    if(getOption("PCMBase.PrintSuffix_u", FALSE)) {
      nLatexGreek <- gsub('Sigma_x', 'Sigma_{u}', nLatexGreek, fixed = TRUE)
      nLatexGreek <- gsub('Sigmae_x', 'Sigma_{e,u}', nLatexGreek, fixed = TRUE)
      nLatexGreek <- gsub('Sigmaj_x', 'Sigma_{j,u}', nLatexGreek, fixed = TRUE)
    } else {
      nLatexGreek <- gsub('Sigma_x', 'Sigma_{x}', nLatexGreek, fixed = TRUE)
      nLatexGreek <- gsub('Sigmae_x', 'Sigma_{e,x}', nLatexGreek, fixed = TRUE)
      nLatexGreek <- gsub('Sigmaj_x', 'Sigma_{j,x}', nLatexGreek, fixed = TRUE)
    }

    nLatexGreek <- gsub('Sigmae', 'Sigma_{e}', nLatexGreek, fixed = TRUE)
    nLatexGreek <- gsub('Sigmaj', 'Sigma_{j}', nLatexGreek, fixed = TRUE)
    nLatexGreek <- gsub('mj', '\\vec{\\mu}_{j}', nLatexGreek, fixed = TRUE)
    if(! (n %in% c('regime', 'type')) ) {
      nLatexGreek <- paste0('$',nLatexGreek,'$')
    }
    nLatexGreek
  }))

  addtorow <- list()
  addtorow$pos <- as.list(seq_len(nrow(x) - 1L))
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
