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

#' @export
PCMParentClasses.OU <- function(model) {
  c("GaussianPCM", "PCM")
}

#' @importFrom expm expm
#' @export
PCMCond.OU <- function(tree, model, r=1, metaI=PCMInfo(NULL, tree, model, verbose), verbose=FALSE) {
  H <- as.matrix(model$H[,, r])
  Theta <- model$Theta[, r]
  Sigma <- as.matrix(model$Sigma_x[,,r]) %*% t(as.matrix(model$Sigma_x[,,r]))
  if(!is.null(model$Sigmae_x)) {
    Sigmae <- as.matrix(model$Sigmae_x[,,r]) %*% t(as.matrix(model$Sigmae_x[,,r]))
  } else {
    Sigmae <- NULL
  }

  V <- PCMCondVOU(H, Sigma, Sigmae, threshold.Lambda_ij = metaI$PCMBase.Threshold.Lambda_ij)
  omega <- function(t, edgeIndex, metaI, e_Ht = NULL) {
    if(is.null(e_Ht)) {
      e_Ht <- expm(-t*H)
    }
    I <- diag(nrow(H))
    (I-e_Ht) %*% Theta
  }
  Phi <- function(t, edgeIndex, metaI, e_Ht = NULL) {
    if(is.null(e_Ht)) {
      expm(-t*H)
    } else {
      e_Ht
    }
  }
  list(omega = omega, Phi = Phi, V = V)
}

#' @export
PCMDescribe.OU <- function(model, ...) {
  "Ornstein-Uhlenbeck branching stochastic process model and a non-phylogenetic variance component"
}

#' @export
PCMSpecifyParams.OU <- function(model, ...) {
  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  list(
    X0 = list(default = rep(0, k),
              type = c("gvector", "full"),
              description = "trait vector at the root; global for all model regimes"),
    H = list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
             type = c("matrix", "full"),
             description = "Selection strength matrix"),
    Theta = list(default = array(0, dim = c(k, R), dimnames = list(NULL, regimes)),
                 type = c("vector", "full"),
                 description = "long-term optimum trait values"),
    Sigma_x = list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                 type = c("matrix", "upper.tri.diag", "positive.diag"),
                 description = "Upper triangular Choleski factor of the unit-time variance-covariance matrix of the BM-process"),
    Sigmae_x = list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "upper.tri.diag", "positive.diag"),
                  description = "Upper triangular Choleski factor of the variance-covariance matrix for the non-phylogenetic trait component"))
}

#' @export
PCMDescribe.OU__posdiagH <- function(model, ...) "OU with positive diagonal matrix H."
#' @export
PCMParentClasses.OU__posdiagH <- function(model) c("OU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.OU__posdiagH <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                 type = c("matrix", "diag", "positive.diag"),
                 description = "Positive diagonal selection strength matrix")
  spec
}

#' @export
PCMDescribe.OU__posdiagH__posdiagSigma_x <- function(model, ...) "OU with positive diagonal matrix H and and with diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.OU__posdiagH__posdiagSigma_x <- function(model) c("OU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.OU__posdiagH__posdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                 type = c("matrix", "diag", "positive.diag"),
                 description = "Positive diagonal selection strength matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec
}

#' @export
PCMDescribe.OU__posdiagH__posdiagSigmae_x <- function(model, ...) "OU with positive diagonal matrix H and and with diagonal Sigmae_x (i.e. no non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.OU__posdiagH__posdiagSigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.OU__posdiagH__posdiagSigmae_x <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                 type = c("matrix", "diag", "positive.diag"),
                 description = "Positive diagonal selection strength matrix")
  spec$Sigmae_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                        type = c("matrix", "diag", "positive.diag"),
                        description = "Positive diagonal Choleski factor of the non-phylogenetic component")
  spec
}

#' @export
PCMDescribe.OU__posdiagH__posdiagSigma_x__posdiagSigmae_x <- function(model, ...) "OU with positive diagonal matrix H and and with positive diagonal Sigma_x and with positive diagonal Sigmae_x (i.e. assuming no phylogenetic or non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.OU__posdiagH__posdiagSigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.OU__posdiagH__posdiagSigmae_x <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                 type = c("matrix", "diag", "positive.diag"),
                 description = "Positive diagonal selection strength matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec$Sigmae_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                        type = c("matrix", "diag", "positive.diag"),
                        description = "Positive diagonal Choleski factor of the non-phylogenetic component")
  spec
}







#' @export
PCMDescribe.OU__noX0 <- function(model, ...) "OU without X0."
#' @export
PCMParentClasses.OU__noX0 <- function(model) c("OU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.OU__noX0 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.OU__noX0__posdiagH <- function(model, ...) "OU without X0 and with positive diagonal matrix H."
#' @export
PCMParentClasses.OU__noX0__posdiagH <- function(model) c("OU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.OU__noX0__posdiagH <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                 type = c("matrix", "diag", "positive.diag"),
                 description = "Positive diagonal selection strength matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.OU__noX0__posdiagH__posdiagSigma_x <- function(model, ...) "OU without X0 and with positive diagonal matrix H and and with diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.OU__noX0__posdiagH__posdiagSigma_x <- function(model) c("OU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.OU__noX0__posdiagH__posdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                 type = c("matrix", "diag", "positive.diag"),
                 description = "Positive diagonal selection strength matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")

  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.OU__noX0__posdiagH__posdiagSigmae_x <- function(model, ...) "OU without X0 and with positive diagonal matrix H and and with diagonal Sigmae_x (i.e. no non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.OU__noX0__posdiagH__posdiagSigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.OU__noX0__posdiagH__posdiagSigmae_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                 type = c("matrix", "diag", "positive.diag"),
                 description = "Positive diagonal selection strength matrix")
  spec$Sigmae_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the non-phylogenetic component")

  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.OU__noX0__posdiagH__posdiagSigma_x__posdiagSigmae_x <- function(model, ...) "OU without X0 and with positive diagonal matrix H and and with positive diagonal Sigma_x and with positive diagonal Sigmae_x (i.e. assuming no phylogenetic or non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.OU__noX0__posdiagH__posdiagSigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.OU__noX0__posdiagH__posdiagSigmae_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                 type = c("matrix", "diag", "positive.diag"),
                 description = "Positive diagonal selection strength matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec$Sigmae_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                        type = c("matrix", "diag", "positive.diag"),
                        description = "Positive diagonal Choleski factor of the non-phylogenetic component")

  spec[!sapply(spec, is.null)]
}




#' @export
PCMDescribe.OU__noSigmae_x <- function(model, ...) "OU without Sigmae_x."
#' @export
PCMParentClasses.OU__noSigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.OU__noSigmae_x <- function(model, ...) {
  spec <- NextMethod()
  spec$Sigmae_x <- NULL
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.OU__noSigmae_x__posdiagH <- function(model, ...) "OU without Sigmae_x and with positive diagonal matrix H."
#' @export
PCMParentClasses.OU__noSigmae_x__posdiagH <- function(model) c("OU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.OU__noSigmae_x__posdiagH <- function(model, ...) {
  spec <- NextMethod()
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                 type = c("matrix", "diag", "positive.diag"),
                 description = "Positive diagonal selection strength matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.OU__noSigmae_x__posdiagH__posdiagSigma_x <- function(model, ...) "OU without Sigmae_x and with positive diagonal matrix H and and with diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.OU__noSigmae_x__posdiagH__posdiagSigma_x <- function(model) c("OU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.OU__noSigmae_x__posdiagH__posdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                 type = c("matrix", "diag", "positive.diag"),
                 description = "Positive diagonal selection strength matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")

  spec[!sapply(spec, is.null)]
}



#' @export
PCMDescribe.OU__noX0__noSigmae_x <- function(model, ...) "OU without X0 and Sigmae_x."
#' @export
PCMParentClasses.OU__noX0__noSigmae_x <- function(model) c("OU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.OU__noX0__noSigmae_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.OU__noX0__noSigmae_x__posdiagH <- function(model, ...) "OU without X0 and Sigmae_x and with positive diagonal matrix H."
#' @export
PCMParentClasses.OU__noX0__noSigmae_x__posdiagH <- function(model) c("OU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.OU__noX0__noSigmae_x__posdiagH <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                 type = c("matrix", "diag", "positive.diag"),
                 description = "Positive diagonal selection strength matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.OU__noX0__noSigmae_x__posdiagH__posdiagSigma_x <- function(model, ...) "OU without X0 and Sigmae_x and with positive diagonal matrix H and and with diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.OU__noX0__noSigmae_x__posdiagH__posdiagSigma_x <- function(model) c("OU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.OU__noX0__noSigmae_x__posdiagH__posdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                 type = c("matrix", "diag", "positive.diag"),
                 description = "Positive diagonal selection strength matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")

  spec[!sapply(spec, is.null)]
}

