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

#' Plot a tree with regimes
#' @importFrom data.table data.table
#' @importFrom grDevices hcl
#' @importFrom ggplot2 aes scale_color_manual
#' @export
PCMTreePlot <- function(tree) {
  if(!require(ggtree)) {
    stop("ERR:02400:PCMBase:Utilities.R:PCMTreePlot:: Calling PCMTreePlot needs ggtree package to be installed from Bioconductor. Check the instructions at https://bioconductor.org/packages/release/bioc/html/ggtree.html. Ggtree was not on CRAN at the time of releasing PCMBase and is not declared as dependency in the PCMBase-description.")
  }

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  if(is.null(tree$edge.regime)) {
    PCMTreeSetDefaultRegime(tree, 1)
  }
  N <- PCMTreeNumTips(tree)
  R <- PCMTreeNumUniqueRegimes(tree)
  palette <- gg_color_hue(R)
  names(palette) <- PCMTreeUniqueRegimes(tree)

  data <- rbind(data.table(node = tree$edge[, 2],
                           regime = as.factor(tree$edge.regime)),
                data.table(node = N+1, regime = as.factor(tree$edge.regime[1])))

  plotTree <- ggtree(tree, layout = 'fan', open.angle = 8, size=.25) %<+% data

  plotTree + aes(color = regime) +
    scale_color_manual(name = "regime", values = palette)
}


#' Scatter plot of 2-dimensional data
#' @importFrom ggplot2 ggplot geom_point scale_size_continuous scale_alpha_continuous geom_text aes theme_gray theme
#' @importFrom data.table data.table is.data.table setkey := setnames
#' @export
PCMPlotTraitData2D <- function(X, tree, labeledTips = NULL) {

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  N <- PCMTreeNumTips(tree)
  R <- PCMTreeNumUniqueRegimes(tree)
  palette <- gg_color_hue(R)
  names(palette) <- PCMTreeUniqueRegimes(tree)

  times <- PCMTreeNodeTimes(tree, tipsOnly = TRUE)
  Xt <- t(X)
  data <- as.data.table(Xt)
  setnames(data, c("x", "y"))
  data[, id:=1:(.N)]
  data[, time:=times]

  data[, regime:=as.factor(sapply(id, function(i) PCMTreeGetRegimeForNode(tree, i)))]
  setkey(data, id)

  if(!is.null(labeledTips)) {
    data[list(labeledTips), label:=id]
  }

  pl <- ggplot(data)
  if(length(unique(PCMTreeNodeTimes(tree, tipsOnly = TRUE))) > 1) {
    # non-ultrametric tree
    pl <- pl + geom_point(aes(x=x, y=y, col=regime, size = time, alpha = time), na.rm = TRUE) +
      scale_size_continuous(range = c(0.2, 2.8)) +
      scale_alpha_continuous(range = c(0.2, 0.75))
  } else {
    pl <- pl + geom_point(aes(x=x, y=y, col=regime), na.rm = TRUE)
  }

  pl <- pl + scale_color_manual(name = "regime", values = palette)

  if(!is.null(labeledTips)) {
    pl <- pl +
      geom_text(aes(label = label, x = x, y = y, col = regime))
  }

  pl
}


#' Extract error information from a formatted error message.
#' @param x character string representing the error message.
#' @description The function searches x for a pattern matching the format
#' 'ERR:5-alphanumeric-character-code:project-name:source-file:error-specifics:'.
#' Specifically it
#' searches for a regular expression pattern "ERR:[0-9a-zA-Z]+:[^:]+:[^:]+:[^:]+:[^:]*:".
#' @return a named list with the parsed error information or NULL, if no match
#' was found. The elements of this list are named as follows:
#' \item{type}{The type of the error message. Usually this is ERROR, but could be
#' WARNING or anything else.}
#' \item{icode}{An an alphanumeric code of the error.}
#' \item{project}{The name of the project locating the code that raised the error.}
#' \item{file}{The name of the source-file containing the code that raised the error.}
#' \item{fun}{The name of the function raising the error}
#' \item{info}{A character string containing additional error-specific information}
#' \item{msg}{A verbal description of the error.}
#' @export
PCMParseErrorMessage <- function(x) {
  res <- try({
    if(is.character(x)) {
      code <- regmatches(x, regexpr("ERR:[0-9a-zA-Z]+:[^:]+:[^:]+:[^:]+:[^:]*:", x))
      if(length(code) > 0) {
        code <- code[1]
        codeL <- strsplit(code, split=":")[[1]]
        list(
          type = codeL[1],
          icode = codeL[2],
          project = codeL[3],
          file = codeL[4],
          fun = codeL[5],
          info = codeL[6],
          code = code,
          msg = x
        )
      } else {
        NULL
      }
    } else {
      NULL
    }
  }, silent = TRUE)

  if(class(res)=="try-error") {
    NULL
  } else {
    res
  }
}

PCMCharacterVectorToRExpression <- function(v) {
  expr <- "c("
  for(i in 1:length(v)) {
    if(i < length(v)) {
      expr <- paste0(expr, "'", v[i], "', ")
    } else {
      # last element
      expr <- paste0(expr, "'", v[i], "')")
    }
  }
  expr
}

#' This function was used to autogenerate the code for PCMSpecify.OU, .JOU, .BM, .White and .DOU
GenerateDefaultParameterizations <- function(libname, pkgname){
  PCMGenerateParameterizations(structure(0.0, class = "OU"), PCMTableParameterizations(structure(0.0, class = "OU"))[1])
  PCMGenerateParameterizations(structure(0.0, class = "DOU"), PCMTableParameterizations(structure(0.0, class = "DOU"))[1])
  PCMGenerateParameterizations(structure(0.0, class = "BM"), PCMTableParameterizations(structure(0.0, class = "BM"))[1])
  PCMGenerateParameterizations(structure(0.0, class = "JOU"), PCMTableParameterizations(structure(0.0, class = "JOU"))[1])
  PCMGenerateParameterizations(structure(0.0, class = "White"), PCMTableParameterizations(structure(0.0, class = "White"))[1])
}

