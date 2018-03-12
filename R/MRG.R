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

#' @name MRG
#' @title Multiple regime Gaussian PCMs

#' @export
is.MRG <- function(x) inherits(x, "MRG")

#' @export
PCMParentClasses.MRG <- function(model) {
  c("GaussianPCM", "PCM")
}

#' @export
PCMDescribe.MRG <- function(model, ...) {
  "Multiple regime Gaussian model"
}

#' @export
PCMCond.MRG <- function(tree, model, r=1, metaI = PCMInfo(NULL, tree, model, verbose), verbose=FALSE) {
  if(! is.MRG(model) ) {
    stop("ERR:02501:PCMBase:MRG.R:PCMCond.MRG:: model should inherit from S3 class 'MRG'.")
  }
  if(!is.character(r)) {
    # integer regime number must be mapped to a character regime because
    # the model object may contain global parameter vectors and matrices apart
    # from models associated with regimes.
    r <- attr(model, "regimes")[r]
  }
  PCMCond(tree, model[[r]], 1, metaI, verbose)
}

