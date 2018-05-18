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
PCMParentClasses.DOU <- function(model) {
  c("GaussianPCM", "PCM")
}

#' @export
PCMCond.DOU <- function(tree, model, r=1, metaI = PCMInfo(NULL, tree, model, verbose), verbose=FALSE) {
  H1 <- as.matrix(model$H1[,,r])
  H2 <- as.matrix(model$H2[,,r])
  Theta <- model$Theta[,r]
  Sigma <- as.matrix(model$Sigma_x[,,r]) %*% t(as.matrix(model$Sigma_x[,,r]))
  if(!is.null(model$Sigmae)) {
    Sigmae <- as.matrix(model$Sigmae_x[,,r]) %*% t(as.matrix(model$Sigmae_x[,,r]))
  } else {
    Sigmae <- NULL
  }

  V <- PCMCondVOU(H2, Sigma, Sigmae, threshold.Lambda_ij = metaI$PCMBase.Threshold.Lambda_ij)
  omega <- function(t, edgeIndex, metaI, e_H1t = NULL) {
    if(is.null(e_H1t)) {
      e_H1t <- expm(-t*H1)
    }
    I <- diag(nrow(H1))
    (I-e_H1t) %*% Theta
  }
  Phi <- function(t, edgeIndex, metaI, e_H1t = NULL) {
    if(is.null(e_H1t)) {
      expm(-t*H1)
    } else {
      e_H1t
    }
  }
  list(omega = omega, Phi = Phi, V = V)
}

#' @export
PCMDescribe.DOU <- function(model, ...) {
  "Double-rate Ornstein-Uhlenbeck branching stochastic process model with a non-phylogenetic variance component"
}

#' @export
PCMSpecifyParams.DOU <- function(model, ...) {
  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  list(
    X0 = list(default = rep(0, k),
              type = c("gvector", "full"),
              description = "trait vector at the root; global for all model regimes"),
    H1 = list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
              type = c("matrix", "full"),
              description = "rate of convergence towards Theta"),
    H2 = list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
              type = c("matrix", "full"),
              description = "rate of decorrelation between pairs of tips"),
    Theta = list(default = array(0, dim = c(k, R), dimnames = list(NULL, regimes)),
                 type = c("vector", "full"),
                 description = "long-term optimum"),
    Sigma_x = list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                 type = c("matrix", "upper.tri.diag", "positive.diag"),
                 description = "Upper triangular Choleski factor of the unit-time variance-covariance matrix of the BM process"),
    Sigmae_x = list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "upper.tri.diag", "positive.diag"),
                  description = "Upper triangular Choleski factor of the variance-covariance matrix for the non-phylogenetic trait component"))
}

#' @export
PCMDescribe.DOU__ZeroH1 <- function(model, ...) "DOU with zero selection strength matrix (H1)."
#' @export
PCMParentClasses.DOU__ZeroH1 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__ZeroH1 <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "fixed"),
                  description = "Fixed zero selection strength matrix")

  spec$Theta = list(default = array(0, dim = c(k, R), dimnames = list(NULL, regimes)),
                    type = c("vector", "fixed"),
                    description = "Fixed 0 long-term optimum (not relevant with zero H1 matrix)")
  spec
}

#' @export
PCMDescribe.DOU__DiagPosdiagH1 <- function(model, ...) "DOU with positive diagonal selection strength matrix (H1)."
#' @export
PCMParentClasses.DOU__DiagPosdiagH1 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__DiagPosdiagH1 <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal rate matrix of convergence towards Theta")
  spec
}


#' @export
PCMDescribe.DOU__DiagPosdiagH2 <- function(model, ...) "DOU with positive diagonal decorrelation rate matrix (H2)."
#' @export
PCMParentClasses.DOU__DiagPosdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__DiagPosdiagH2 <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec
}

#' @export
PCMDescribe.DOU__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) "DOU with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec
}

#' @export
PCMDescribe.DOU__DiagPosdiagH2__DiagPosdiagSigmae_x <- function(model, ...) "DOU with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigmae_x (i.e. no non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__DiagPosdiagH2__DiagPosdiagSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__DiagPosdiagH2__DiagPosdiagSigmae_x <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec$Sigmae_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the non-phylogenetic variance-covariance matrix")
  spec
}

#' @export
PCMDescribe.DOU__DiagPosdiagH2__DiagPosdiagSigma_x__DiagPosdiagSigmae_x <- function(model, ...) "DOU with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x and Sigmae_x (i.e. no phylogenetic or non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__DiagPosdiagH2__DiagPosdiagSigma_x__DiagPosdiagSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__DiagPosdiagH2__DiagPosdiagSigma_x__DiagPosdiagSigmae_x <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec$Sigmae_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                        type = c("matrix", "diag", "positive.diag"),
                        description = "Positive diagonal Choleski factor of the non-phylogenetic variance-covariance matrix")
  spec
}







#' @export
PCMDescribe.DOU__ZeroH1__DiagPosdiagH2 <- function(model, ...) "DOU with zero selection strength matrix (H1) and positive diagonal decorrelation rate matrix (H2)."
#' @export
PCMParentClasses.DOU__ZeroH1__DiagPosdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__ZeroH1__DiagPosdiagH2 <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "fixed"),
                  description = "Fixed zero selection strength matrix")

  spec$Theta = list(default = array(0, dim = c(k, R), dimnames = list(NULL, regimes)),
                    type = c("vector", "fixed"),
                    description = "Fixed 0 long-term optimum (not relevant with zero H1 matrix)")
  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec
}

#' @export
PCMDescribe.DOU__ZeroH1__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) "DOU with zero selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__ZeroH1__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__ZeroH1__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "fixed"),
                  description = "Fixed zero selection strength matrix")

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec
}

#' @export
PCMDescribe.DOU__ZeroH1__DiagPosdiagH2__DiagPosdiagSigmae_x <- function(model, ...) "DOU with zero selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigmae_x (i.e. no non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__ZeroH1__DiagPosdiagH2__DiagPosdiagSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__ZeroH1__DiagPosdiagH2__DiagPosdiagSigmae_x <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "fixed"),
                  description = "Fixed zero selection strength matrix")

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec$Sigmae_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                        type = c("matrix", "diag", "positive.diag"),
                        description = "Positive diagonal Choleski factor of the non-phylogenetic variance-covariance matrix")
  spec
}

#' @export
PCMDescribe.DOU__ZeroH1__DiagPosdiagH2__DiagPosdiagSigma_x__DiagPosdiagSigmae_x <- function(model, ...) "DOU with zero selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x and Sigmae_x (i.e. no phylogenetic or non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__ZeroH1__DiagPosdiagH2__DiagPosdiagSigma_x__DiagPosdiagSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__ZeroH1__DiagPosdiagH2__DiagPosdiagSigma_x__DiagPosdiagSigmae_x <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "fixed"),
                  description = "Fixed zero selection strength matrix")

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec$Sigmae_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                        type = c("matrix", "diag", "positive.diag"),
                        description = "Positive diagonal Choleski factor of the non-phylogenetic variance-covariance matrix")
  spec
}













#' @export
PCMDescribe.DOU__DiagPosdiagH1__DiagPosdiagH2 <- function(model, ...) "DOU with positive diagonal selection strength matrix (H1) and positive diagonal decorrelation rate matrix (H2)."
#' @export
PCMParentClasses.DOU__DiagPosdiagH1__DiagPosdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__DiagPosdiagH1__DiagPosdiagH2 <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Fixed zero selection strength matrix")

  spec$Theta = list(default = array(0, dim = c(k, R), dimnames = list(NULL, regimes)),
                    type = c("vector", "fixed"),
                    description = "Fixed 0 long-term optimum (not relevant with zero H1 matrix)")
  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec
}

#' @export
PCMDescribe.DOU__DiagPosdiagH1__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) "DOU with positive diagonal selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__DiagPosdiagH1__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__DiagPosdiagH1__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Fixed zero selection strength matrix")

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec
}

#' @export
PCMDescribe.DOU__DiagPosdiagH1__DiagPosdiagH2__DiagPosdiagSigmae_x <- function(model, ...) "DOU with positive diagonal selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigmae_x (i.e. no non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__DiagPosdiagH1__DiagPosdiagH2__DiagPosdiagSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__DiagPosdiagH1__DiagPosdiagH2__DiagPosdiagSigmae_x <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Fixed zero selection strength matrix")

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec$Sigmae_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                        type = c("matrix", "diag", "positive.diag"),
                        description = "Positive diagonal Choleski factor of the non-phylogenetic variance-covariance matrix")
  spec
}

#' @export
PCMDescribe.DOU__DiagPosdiagH1__DiagPosdiagH2__DiagPosdiagSigma_x__DiagPosdiagSigmae_x <- function(model, ...) "DOU with positive diagonal selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x and Sigmae_x (i.e. no phylogenetic or non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__DiagPosdiagH1__DiagPosdiagH2__DiagPosdiagSigma_x__DiagPosdiagSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__DiagPosdiagH1__DiagPosdiagH2__DiagPosdiagSigma_x__DiagPosdiagSigmae_x <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Fixed zero selection strength matrix")

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec$Sigmae_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                        type = c("matrix", "diag", "positive.diag"),
                        description = "Positive diagonal Choleski factor of the non-phylogenetic variance-covariance matrix")
  spec
}



















##########################
### noX0



#' @export
PCMDescribe.DOU__NoX0 <- function(model, ...) "DOU without X0."
#' @export
PCMParentClasses.DOU__NoX0 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec[!sapply(spec, is.null)]
}



#' @export
PCMDescribe.DOU__NoX0__ZeroH1 <- function(model, ...) "DOU without X0 and with zero selection strength matrix (H1)."
#' @export
PCMParentClasses.DOU__NoX0__ZeroH1 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__ZeroH1 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "fixed"),
                  description = "Fixed zero selection strength matrix")

  spec$Theta = list(default = array(0, dim = c(k, R), dimnames = list(NULL, regimes)),
                    type = c("vector", "fixed"),
                    description = "Fixed 0 long-term optimum (not relevant with zero H1 matrix)")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoX0__DiagPosdiagH1 <- function(model, ...) "DOU without X0 and with positive diagonal selection strength matrix (H1)."
#' @export
PCMParentClasses.DOU__NoX0__DiagPosdiagH1 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__DiagPosdiagH1 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal rate matrix of convergence towards Theta")
  spec[!sapply(spec, is.null)]
}


#' @export
PCMDescribe.DOU__NoX0__DiagPosdiagH2 <- function(model, ...) "DOU without X0 and with positive diagonal decorrelation rate matrix (H2)."
#' @export
PCMParentClasses.DOU__NoX0__DiagPosdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__DiagPosdiagH2 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoX0__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) "DOU without X0 and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__NoX0__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoX0__DiagPosdiagH2__DiagPosdiagSigmae_x <- function(model, ...) "DOU without X0 and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigmae_x (i.e. no non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__NoX0__DiagPosdiagH2__DiagPosdiagSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__DiagPosdiagH2__DiagPosdiagSigmae_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec$Sigmae_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                        type = c("matrix", "diag", "positive.diag"),
                        description = "Positive diagonal Choleski factor of the non-phylogenetic variance-covariance matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoX0__DiagPosdiagH2__DiagPosdiagSigma_x__DiagPosdiagSigmae_x <- function(model, ...) "DOU without X0 and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x and Sigmae_x (i.e. no phylogenetic or non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__NoX0__DiagPosdiagH2__DiagPosdiagSigma_x__DiagPosdiagSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__DiagPosdiagH2__DiagPosdiagSigma_x__DiagPosdiagSigmae_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec$Sigmae_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                        type = c("matrix", "diag", "positive.diag"),
                        description = "Positive diagonal Choleski factor of the non-phylogenetic variance-covariance matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoX0__ZeroH1__DiagPosdiagH2 <- function(model, ...) "DOU without X0 and with zero selection strength matrix (H1) and positive diagonal decorrelation rate matrix (H2)."
#' @export
PCMParentClasses.DOU__NoX0__ZeroH1__DiagPosdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__ZeroH1__DiagPosdiagH2 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "fixed"),
                  description = "Fixed zero selection strength matrix")

  spec$Theta = list(default = array(0, dim = c(k, R), dimnames = list(NULL, regimes)),
                    type = c("vector", "fixed"),
                    description = "Fixed 0 long-term optimum (not relevant with zero H1 matrix)")
  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoX0__ZeroH1__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) "DOU without X0 and with zero selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__NoX0__ZeroH1__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__ZeroH1__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "fixed"),
                  description = "Fixed zero selection strength matrix")

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoX0__ZeroH1__DiagPosdiagH2__DiagPosdiagSigmae_x <- function(model, ...) "DOU without X0 and with zero selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigmae_x (i.e. no non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__NoX0__ZeroH1__DiagPosdiagH2__DiagPosdiagSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__ZeroH1__DiagPosdiagH2__DiagPosdiagSigmae_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "fixed"),
                  description = "Fixed zero selection strength matrix")

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec$Sigmae_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                        type = c("matrix", "diag", "positive.diag"),
                        description = "Positive diagonal Choleski factor of the non-phylogenetic variance-covariance matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoX0__ZeroH1__DiagPosdiagH2__DiagPosdiagSigma_x__DiagPosdiagSigmae_x <- function(model, ...) "DOU without X0 and with zero selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x and Sigmae_x (i.e. no phylogenetic or non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__NoX0__ZeroH1__DiagPosdiagH2__DiagPosdiagSigma_x__DiagPosdiagSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__ZeroH1__DiagPosdiagH2__DiagPosdiagSigma_x__DiagPosdiagSigmae_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "fixed"),
                  description = "Fixed zero selection strength matrix")

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec$Sigmae_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                        type = c("matrix", "diag", "positive.diag"),
                        description = "Positive diagonal Choleski factor of the non-phylogenetic variance-covariance matrix")
  spec[!sapply(spec, is.null)]
}


#' @export
PCMDescribe.DOU__NoX0__DiagPosdiagH1__DiagPosdiagH2 <- function(model, ...) "DOU without X0 and with positive diagonal selection strength matrix (H1) and positive diagonal decorrelation rate matrix (H2)."
#' @export
PCMParentClasses.DOU__NoX0__DiagPosdiagH1__DiagPosdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__DiagPosdiagH1__DiagPosdiagH2 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Fixed zero selection strength matrix")

  spec$Theta = list(default = array(0, dim = c(k, R), dimnames = list(NULL, regimes)),
                    type = c("vector", "fixed"),
                    description = "Fixed 0 long-term optimum (not relevant with zero H1 matrix)")
  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoX0__DiagPosdiagH1__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) "DOU without X0 and with positive diagonal selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__NoX0__DiagPosdiagH1__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__DiagPosdiagH1__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Fixed zero selection strength matrix")

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoX0__DiagPosdiagH1__DiagPosdiagH2__DiagPosdiagSigmae_x <- function(model, ...) "DOU without X0 and with positive diagonal selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigmae_x (i.e. no non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__NoX0__DiagPosdiagH1__DiagPosdiagH2__DiagPosdiagSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__DiagPosdiagH1__DiagPosdiagH2__DiagPosdiagSigmae_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Fixed zero selection strength matrix")

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec$Sigmae_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                        type = c("matrix", "diag", "positive.diag"),
                        description = "Positive diagonal Choleski factor of the non-phylogenetic variance-covariance matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoX0__DiagPosdiagH1__DiagPosdiagH2__DiagPosdiagSigma_x__DiagPosdiagSigmae_x <- function(model, ...) "DOU without X0 and with positive diagonal selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x and Sigmae_x (i.e. no phylogenetic or non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__NoX0__DiagPosdiagH1__DiagPosdiagH2__DiagPosdiagSigma_x__DiagPosdiagSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__DiagPosdiagH1__DiagPosdiagH2__DiagPosdiagSigma_x__DiagPosdiagSigmae_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Fixed zero selection strength matrix")

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec$Sigmae_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                        type = c("matrix", "diag", "positive.diag"),
                        description = "Positive diagonal Choleski factor of the non-phylogenetic variance-covariance matrix")
  spec[!sapply(spec, is.null)]
}
























#####################
#### no Sigmae_x


#' @export
PCMDescribe.DOU__NoSigmae_x <- function(model, ...) "DOU without Sigmae_x."
#' @export
PCMParentClasses.DOU__NoSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoSigmae_x <- function(model, ...) {
  spec <- NextMethod()
  spec$Sigmae_x <- NULL
  spec[!sapply(spec, is.null)]
}


#' @export
PCMDescribe.DOU__NoSigmae_x__ZeroH1 <- function(model, ...) "DOU without Sigmae_x and with zero selection strength matrix (H1)."
#' @export
PCMParentClasses.DOU__NoSigmae_x__ZeroH1 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoSigmae_x__ZeroH1 <- function(model, ...) {
  spec <- NextMethod()
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "fixed"),
                  description = "Fixed zero selection strength matrix")

  spec$Theta = list(default = array(0, dim = c(k, R), dimnames = list(NULL, regimes)),
                    type = c("vector", "fixed"),
                    description = "Fixed 0 long-term optimum (not relevant with zero H1 matrix)")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoSigmae_x__DiagPosdiagH1 <- function(model, ...) "DOU without Sigmae_x and with positive diagonal selection strength matrix (H1)."
#' @export
PCMParentClasses.DOU__NoSigmae_x__DiagPosdiagH1 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoSigmae_x__DiagPosdiagH1 <- function(model, ...) {
  spec <- NextMethod()
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal rate matrix of convergence towards Theta")
  spec[!sapply(spec, is.null)]
}


#' @export
PCMDescribe.DOU__NoSigmae_x__DiagPosdiagH2 <- function(model, ...) "DOU without Sigmae_x and with positive diagonal decorrelation rate matrix (H2)."
#' @export
PCMParentClasses.DOU__NoSigmae_x__DiagPosdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoSigmae_x__DiagPosdiagH2 <- function(model, ...) {
  spec <- NextMethod()
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoSigmae_x__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) "DOU without Sigmae_x and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__NoSigmae_x__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoSigmae_x__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec[!sapply(spec, is.null)]
}


#' @export
PCMDescribe.DOU__NoSigmae_x__ZeroH1__DiagPosdiagH2 <- function(model, ...) "DOU without Sigmae_x and with zero selection strength matrix (H1) and positive diagonal decorrelation rate matrix (H2)."
#' @export
PCMParentClasses.DOU__NoSigmae_x__ZeroH1__DiagPosdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoSigmae_x__ZeroH1__DiagPosdiagH2 <- function(model, ...) {
  spec <- NextMethod()
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "fixed"),
                  description = "Fixed zero selection strength matrix")

  spec$Theta = list(default = array(0, dim = c(k, R), dimnames = list(NULL, regimes)),
                    type = c("vector", "fixed"),
                    description = "Fixed 0 long-term optimum (not relevant with zero H1 matrix)")
  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoSigmae_x__ZeroH1__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) "DOU without Sigmae_x and with zero selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__NoSigmae_x__ZeroH1__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoSigmae_x__ZeroH1__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "fixed"),
                  description = "Fixed zero selection strength matrix")

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoSigmae_x__DiagPosdiagH1__DiagPosdiagH2 <- function(model, ...) "DOU without Sigmae_x and with positive diagonal selection strength matrix (H1) and positive diagonal decorrelation rate matrix (H2)."
#' @export
PCMParentClasses.DOU__NoSigmae_x__DiagPosdiagH1__DiagPosdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoSigmae_x__DiagPosdiagH1__DiagPosdiagH2 <- function(model, ...) {
  spec <- NextMethod()
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Fixed zero selection strength matrix")

  spec$Theta = list(default = array(0, dim = c(k, R), dimnames = list(NULL, regimes)),
                    type = c("vector", "fixed"),
                    description = "Fixed 0 long-term optimum (not relevant with zero H1 matrix)")
  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoSigmae_x__DiagPosdiagH1__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) "DOU without Sigmae_x and with positive diagonal selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__NoSigmae_x__DiagPosdiagH1__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoSigmae_x__DiagPosdiagH1__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Fixed zero selection strength matrix")

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec[!sapply(spec, is.null)]
}




#######################
### noX0__NoSigmae_x

#' @export
PCMDescribe.DOU__NoX0__NoSigmae_x <- function(model, ...) "DOU without X0 and Sigmae_x."
#' @export
PCMParentClasses.DOU__NoX0__NoSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__NoSigmae_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL
  spec[!sapply(spec, is.null)]
}


#' @export
PCMDescribe.DOU__NoX0__NoSigmae_x__ZeroH1 <- function(model, ...) "DOU without X0 and without Sigmae_x and with zero selection strength matrix (H1)."
#' @export
PCMParentClasses.DOU__NoX0__NoSigmae_x__ZeroH1 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__NoSigmae_x__ZeroH1 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "fixed"),
                  description = "Fixed zero selection strength matrix")

  spec$Theta = list(default = array(0, dim = c(k, R), dimnames = list(NULL, regimes)),
                    type = c("vector", "fixed"),
                    description = "Fixed 0 long-term optimum (not relevant with zero H1 matrix)")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoX0__NoSigmae_x__PosdiagH1 <- function(model, ...) "DOU without X0 and without Sigmae_x and with positive diagonal elements of the (non-diagonal) selection strength matrix (H1)."
#' @export
PCMParentClasses.DOU__NoX0__NoSigmae_x__PosdiagH1 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__NoSigmae_x__PosdiagH1 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "full", "positive.diag"),
                  description = "Rate matrix of convergence towards Theta with positive diagonal elements")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoX0__NoSigmae_x__SymmetricPosdiagH1 <- function(model, ...) "DOU without X0 and without Sigmae_x and with symmetric selection strength matrix (H1) with positive diagonal elements."
#' @export
PCMParentClasses.DOU__NoX0__NoSigmae_x__SymmetricPosdiagH1 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__NoSigmae_x__SymmetricPosdiagH1 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "symmetric", "positive.diag"),
                  description = "Symmetric rate matrix of convergence towards Theta with positive diagonal elements")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoX0__NoSigmae_x__DiagPosdiagH1 <- function(model, ...) "DOU without X0 and without Sigmae_x and with positive diagonal selection strength matrix (H1)."
#' @export
PCMParentClasses.DOU__NoX0__NoSigmae_x__DiagPosdiagH1 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__NoSigmae_x__DiagPosdiagH1 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal rate matrix of convergence towards Theta")
  spec[!sapply(spec, is.null)]
}


#' @export
PCMDescribe.DOU__NoX0__NoSigmae_x__PosdiagH2 <- function(model, ...) "DOU without X0 and without Sigmae_x and with positive diagonal elements of the (non-diagonal) decorrelation rate matrix (H2)."
#' @export
PCMParentClasses.DOU__NoX0__NoSigmae_x__PosdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__NoSigmae_x__PosdiagH2 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "full", "positive.diag"),
                  description = "Decorrelation rate matrix with positive diagonal elements")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoX0__NoSigmae_x__SymmetricPosdiagH2 <- function(model, ...) "DOU without X0 and without Sigmae_x and with symmetric decorrelation rate matrix (H2) with positive diagonal elements."
#' @export
PCMParentClasses.DOU__NoX0__NoSigmae_x__SymmetricPosdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__NoSigmae_x__SymmetricPosdiagH2 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "symmetric", "positive.diag"),
                  description = "Symmetric decorrelation rate matrix with positive diagonal elements")
  spec[!sapply(spec, is.null)]
}


#' @export
PCMDescribe.DOU__NoX0__NoSigmae_x__DiagPosdiagH2 <- function(model, ...) "DOU without X0 and without Sigmae_x and with positive diagonal decorrelation rate matrix (H2)."
#' @export
PCMParentClasses.DOU__NoX0__NoSigmae_x__DiagPosdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__NoSigmae_x__DiagPosdiagH2 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoX0__NoSigmae_x__PosdiagH2__DiagPosdiagSigma_x <- function(model, ...) "DOU without X0 and without Sigmae_x and with positive diagonal elements of the (non-diagonal) decorrelation rate matrix (H2) and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__NoX0__NoSigmae_x__PosdiagH2__DiagPosdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__NoSigmae_x__PosdiagH2__DiagPosdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "full", "positive.diag"),
                  description = "Decorrelation rate matrix with positive diagonal elements")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoX0__NoSigmae_x__SymmetricPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) "DOU without X0 and without Sigmae_x and with symmetric decorrelation rate matrix (H2) with positive diagonal elements and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__NoX0__NoSigmae_x__SymmetricPosdiagH2__DiagPosdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__NoSigmae_x__SymmetricPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "symmetric", "positive.diag"),
                  description = "Symmetric decorrelation rate matrix with positive diagonal elements")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoX0__NoSigmae_x__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) "DOU without X0 and without Sigmae_x and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__NoX0__NoSigmae_x__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__NoSigmae_x__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec[!sapply(spec, is.null)]
}


#' @export
PCMDescribe.DOU__NoX0__NoSigmae_x__ZeroH1__PosdiagH2 <- function(model, ...) "DOU without X0 and without Sigmae_x and with zero selection strength matrix (H1) and decorrelation rate matrix (H2) with positive diagonal elements."
#' @export
PCMParentClasses.DOU__NoX0__NoSigmae_x__ZeroH1__PosdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__NoSigmae_x__ZeroH1__PosdiagH2 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "fixed"),
                  description = "Fixed zero selection strength matrix")

  spec$Theta = list(default = array(0, dim = c(k, R), dimnames = list(NULL, regimes)),
                    type = c("vector", "fixed"),
                    description = "Fixed 0 long-term optimum (not relevant with zero H1 matrix)")
  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "full", "positive.diag"),
                  description = "Decorrelation rate matrix with positive diagonal elements")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoX0__NoSigmae_x__ZeroH1__SymmetricPosdiagH2 <- function(model, ...) "DOU without X0 and without Sigmae_x and with zero selection strength matrix (H1) and symmetric decorrelation rate matrix (H2) with positive diagonal elements."
#' @export
PCMParentClasses.DOU__NoX0__NoSigmae_x__ZeroH1__SymmetricPosdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__NoSigmae_x__ZeroH1__SymmetricPosdiagH2 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "fixed"),
                  description = "Fixed zero selection strength matrix")

  spec$Theta = list(default = array(0, dim = c(k, R), dimnames = list(NULL, regimes)),
                    type = c("vector", "fixed"),
                    description = "Fixed 0 long-term optimum (not relevant with zero H1 matrix)")
  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "symmetric", "positive.diag"),
                  description = "Symmetric decorrelation rate matrix with positive diagonal elements")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoX0__NoSigmae_x__ZeroH1__DiagPosdiagH2 <- function(model, ...) "DOU without X0 and without Sigmae_x and with zero selection strength matrix (H1) and positive diagonal decorrelation rate matrix (H2)."
#' @export
PCMParentClasses.DOU__NoX0__NoSigmae_x__ZeroH1__DiagPosdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__NoSigmae_x__ZeroH1__DiagPosdiagH2 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "fixed"),
                  description = "Fixed zero selection strength matrix")

  spec$Theta = list(default = array(0, dim = c(k, R), dimnames = list(NULL, regimes)),
                    type = c("vector", "fixed"),
                    description = "Fixed 0 long-term optimum (not relevant with zero H1 matrix)")
  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoX0__NoSigmae_x__ZeroH1__PosdiagH2__DiagPosdiagSigma_x <- function(model, ...) "DOU without X0 and without Sigmae_x and with zero selection strength matrix (H1) and with a decorrelation rate matrix (H2) with positive diagonal elements and with positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__NoX0__NoSigmae_x__ZeroH1__PosdiagH2__DiagPosdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__NoSigmae_x__ZeroH1__PosdiagH2__DiagPosdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "fixed"),
                  description = "Fixed zero selection strength matrix")

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "full", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoX0__NoSigmae_x__ZeroH1__SymmetricPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) "DOU without X0 and without Sigmae_x and with zero selection strength matrix (H1) and with a symmetric decorrelation rate matrix (H2) with positive diagonal elements and with positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__NoX0__NoSigmae_x__ZeroH1__SymmetricPosdiagH2__DiagPosdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__NoSigmae_x__ZeroH1__SymmetricPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "fixed"),
                  description = "Fixed zero selection strength matrix")

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "symmetric", "positive.diag"),
                  description = "Symmetric decorrelation rate matrix with positive diagonal elements")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec[!sapply(spec, is.null)]
}



#' @export
PCMDescribe.DOU__NoX0__NoSigmae_x__ZeroH1__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) "DOU without X0 and without Sigmae_x and with zero selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__NoX0__NoSigmae_x__ZeroH1__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__NoSigmae_x__ZeroH1__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "fixed"),
                  description = "Fixed zero selection strength matrix")

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec[!sapply(spec, is.null)]
}


#' @export
PCMDescribe.DOU__NoX0__NoSigmae_x__DiagPosdiagH1__PosdiagH2 <- function(model, ...) "DOU without X0 and without Sigmae_x and with positive diagonal selection strength matrix (H1) and decorrelation rate matrix (H2) with positive diagonal elements."
#' @export
PCMParentClasses.DOU__NoX0__NoSigmae_x__DiagPosdiagH1__PosdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__NoSigmae_x__DiagPosdiagH1__PosdiagH2 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Fixed zero selection strength matrix")

  spec$Theta = list(default = array(0, dim = c(k, R), dimnames = list(NULL, regimes)),
                    type = c("vector", "fixed"),
                    description = "Fixed 0 long-term optimum (not relevant with zero H1 matrix)")
  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "full", "positive.diag"),
                  description = "Decorrelation rate matrix with positive diagonal elements")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoX0__NoSigmae_x__DiagPosdiagH1__SymmetricPosdiagH2 <- function(model, ...) "DOU without X0 and without Sigmae_x and with positive diagonal selection strength matrix (H1) and symmetric decorrelation rate matrix (H2) with positive diagonal elements."
#' @export
PCMParentClasses.DOU__NoX0__NoSigmae_x__DiagPosdiagH1__SymmetricPosdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__NoSigmae_x__DiagPosdiagH1__SymmetricPosdiagH2 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Fixed zero selection strength matrix")

  spec$Theta = list(default = array(0, dim = c(k, R), dimnames = list(NULL, regimes)),
                    type = c("vector", "fixed"),
                    description = "Fixed 0 long-term optimum (not relevant with zero H1 matrix)")
  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "symmetric", "positive.diag"),
                  description = "Symmetric decorrelation rate matrix with positive diagonal elements")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoX0__NoSigmae_x__DiagPosdiagH1__DiagPosdiagH2 <- function(model, ...) "DOU without X0 and without Sigmae_x and with positive diagonal selection strength matrix (H1) and positive diagonal decorrelation rate matrix (H2)."
#' @export
PCMParentClasses.DOU__NoX0__NoSigmae_x__DiagPosdiagH1__DiagPosdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__NoSigmae_x__DiagPosdiagH1__DiagPosdiagH2 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Fixed zero selection strength matrix")

  spec$Theta = list(default = array(0, dim = c(k, R), dimnames = list(NULL, regimes)),
                    type = c("vector", "fixed"),
                    description = "Fixed 0 long-term optimum (not relevant with zero H1 matrix)")
  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoX0__NoSigmae_x__DiagPosdiagH1__PosdiagH2__DiagPosdiagSigma_x <- function(model, ...) "DOU without X0 and without Sigmae_x and with positive diagonal selection strength matrix (H1) and with decorrelation rate matrix (H2) with positive diagonal elements and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__NoX0__NoSigmae_x__DiagPosdiagH1__PosdiagH2__DiagPosdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__NoSigmae_x__DiagPosdiagH1__PosdiagH2__DiagPosdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Fixed zero selection strength matrix")

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "full", "positive.diag"),
                  description = "Decorrelation rate matrix with positive diagonal elements")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec[!sapply(spec, is.null)]
}


#' @export
PCMDescribe.DOU__NoX0__NoSigmae_x__DiagPosdiagH1__SymmetricPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) "DOU without X0 and without Sigmae_x and with positive diagonal selection strength matrix (H1) and with symmetric decorrelation rate matrix (H2) with positive diagonal elements and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__NoX0__NoSigmae_x__DiagPosdiagH1__SymmetricPosdiagH2__DiagPosdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__NoSigmae_x__DiagPosdiagH1__SymmetricPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Fixed zero selection strength matrix")

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "symmetric", "positive.diag"),
                  description = "Symmetric decorrelation rate matrix with positive diagonal elements")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__NoX0__NoSigmae_x__DiagPosdiagH1__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) "DOU without X0 and without Sigmae_x and with positive diagonal selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__NoX0__NoSigmae_x__DiagPosdiagH1__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__NoX0__NoSigmae_x__DiagPosdiagH1__DiagPosdiagH2__DiagPosdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Fixed zero selection strength matrix")

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Positive diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec[!sapply(spec, is.null)]
}
