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

#' @name GG
#' @title Generalized Gaussian PCMs

PCMDescribe.GG <- function(model, ...) {
  R <- length(attr(model, "regimes", exact = TRUE))
  paste0("Generalized Gaussian model with ", PCMN)
}
#' @describeIn GG
#' @inheritParams PCMCond
#' @export
PCMCond.GG <- function(tree, model, r=1, metaI = PCMInfo(NULL, tree, model, verbose), verbose=FALSE) {
  if(! is.GG(model) ) {
    stop("ERR:02501:PCMBase:GG.R:PCMCond.GG:: model should inherit from S3 class 'GG'.")
  }
  if(!is.character(r)) {
    # integer regime number must be mapped to a character regime because
    # the model object may contain global parameter vectors and matrices apart
    # from models associated with regimes.
    r <- attr(model, "regimes")[r]
  }
  PCMCond(tree, model[[r]], 1, metaI, verbose)
}

#' @describeIn GG
#' @inheritParams is.PCM
#'
#' @export
is.GG <- function(x) inherits(x, "GG")

#' @describeIn GG
#' @inheritParams format.PCM
#' @export
format.GG <- function(x, ...) {
  res <- format.PCM(x, ...)
  res <- paste0(res, "Regimes:\n")
  for(i in 1:length(x)) {
    if(is.PCM(x[[i]])) {
      if(!is.null(names(x))) {
        res <- paste0(res, names(x)[i], ": ")
      }
      res <- paste0(res, format(x[[i]], ...), "\n")
    }
  }
}

