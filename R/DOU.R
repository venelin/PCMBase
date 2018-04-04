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
PCMDescribe.DOU__zeroH1 <- function(model, ...) "DOU with zero selection strength matrix (H1)."
#' @export
PCMParentClasses.DOU__zeroH1 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__zeroH1 <- function(model, ...) {
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
PCMDescribe.DOU__posdiagH1 <- function(model, ...) "DOU with positive diagonal selection strength matrix (H1)."
#' @export
PCMParentClasses.DOU__posdiagH1 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__posdiagH1 <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec
}


#' @export
PCMDescribe.DOU__posdiagH2 <- function(model, ...) "DOU with positive diagonal decorrelation rate matrix (H2)."
#' @export
PCMParentClasses.DOU__posdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__posdiagH2 <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec
}

#' @export
PCMDescribe.DOU__posdiagH2__posdiagSigma_x <- function(model, ...) "DOU with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__posdiagH2__posdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__posdiagH2__posdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec
}

#' @export
PCMDescribe.DOU__posdiagH2__posdiagSigmae_x <- function(model, ...) "DOU with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigmae_x (i.e. no non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__posdiagH2__posdiagSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__posdiagH2__posdiagSigmae_x <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec$Sigmae_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the non-phylogenetic variance-covariance matrix")
  spec
}

#' @export
PCMDescribe.DOU__posdiagH2__posdiagSigma_x__posdiagSigmae_x <- function(model, ...) "DOU with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x and Sigmae_x (i.e. no phylogenetic or non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__posdiagH2__posdiagSigma_x__posdiagSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__posdiagH2__posdiagSigma_x__posdiagSigmae_x <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec$Sigmae_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                        type = c("matrix", "diag", "positive.diag"),
                        description = "Positive diagonal Choleski factor of the non-phylogenetic variance-covariance matrix")
  spec
}







#' @export
PCMDescribe.DOU__zeroH1__posdiagH2 <- function(model, ...) "DOU with zero selection strength matrix (H1) and positive diagonal decorrelation rate matrix (H2)."
#' @export
PCMParentClasses.DOU__zeroH1__posdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__zeroH1__posdiagH2 <- function(model, ...) {
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
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec
}

#' @export
PCMDescribe.DOU__zeroH1__posdiagH2__posdiagSigma_x <- function(model, ...) "DOU with zero selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__zeroH1__posdiagH2__posdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__zeroH1__posdiagH2__posdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "fixed"),
                  description = "Fixed zero selection strength matrix")

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec
}

#' @export
PCMDescribe.DOU__zeroH1__posdiagH2__posdiagSigmae_x <- function(model, ...) "DOU with zero selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigmae_x (i.e. no non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__zeroH1__posdiagH2__posdiagSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__zeroH1__posdiagH2__posdiagSigmae_x <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "fixed"),
                  description = "Fixed zero selection strength matrix")

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec$Sigmae_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                        type = c("matrix", "diag", "positive.diag"),
                        description = "Positive diagonal Choleski factor of the non-phylogenetic variance-covariance matrix")
  spec
}

#' @export
PCMDescribe.DOU__zeroH1__posdiagH2__posdiagSigma_x__posdiagSigmae_x <- function(model, ...) "DOU with zero selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x and Sigmae_x (i.e. no phylogenetic or non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__zeroH1__posdiagH2__posdiagSigma_x__posdiagSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__zeroH1__posdiagH2__posdiagSigma_x__posdiagSigmae_x <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "fixed"),
                  description = "Fixed zero selection strength matrix")

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec$Sigmae_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                        type = c("matrix", "diag", "positive.diag"),
                        description = "Positive diagonal Choleski factor of the non-phylogenetic variance-covariance matrix")
  spec
}













#' @export
PCMDescribe.DOU__posdiagH1__posdiagH2 <- function(model, ...) "DOU with positive diagonal selection strength matrix (H1) and positive diagonal decorrelation rate matrix (H2)."
#' @export
PCMParentClasses.DOU__posdiagH1__posdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__posdiagH1__posdiagH2 <- function(model, ...) {
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
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec
}

#' @export
PCMDescribe.DOU__posdiagH1__posdiagH2__posdiagSigma_x <- function(model, ...) "DOU with positive diagonal selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__posdiagH1__posdiagH2__posdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__posdiagH1__posdiagH2__posdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Fixed zero selection strength matrix")

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec
}

#' @export
PCMDescribe.DOU__posdiagH1__posdiagH2__posdiagSigmae_x <- function(model, ...) "DOU with positive diagonal selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigmae_x (i.e. no non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__posdiagH1__posdiagH2__posdiagSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__posdiagH1__posdiagH2__posdiagSigmae_x <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Fixed zero selection strength matrix")

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec$Sigmae_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                        type = c("matrix", "diag", "positive.diag"),
                        description = "Positive diagonal Choleski factor of the non-phylogenetic variance-covariance matrix")
  spec
}

#' @export
PCMDescribe.DOU__posdiagH1__posdiagH2__posdiagSigma_x__posdiagSigmae_x <- function(model, ...) "DOU with positive diagonal selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x and Sigmae_x (i.e. no phylogenetic or non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__posdiagH1__posdiagH2__posdiagSigma_x__posdiagSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__posdiagH1__posdiagH2__posdiagSigma_x__posdiagSigmae_x <- function(model, ...) {
  spec <- NextMethod()

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Fixed zero selection strength matrix")

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Non-negative diagonal decorrelation rate matrix")
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
PCMDescribe.DOU__noX0 <- function(model, ...) "DOU without X0."
#' @export
PCMParentClasses.DOU__noX0 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noX0 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec[!sapply(spec, is.null)]
}



#' @export
PCMDescribe.DOU__noX0__zeroH1 <- function(model, ...) "DOU without X0 and with zero selection strength matrix (H1)."
#' @export
PCMParentClasses.DOU__noX0__zeroH1 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noX0__zeroH1 <- function(model, ...) {
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
PCMDescribe.DOU__noX0__posdiagH1 <- function(model, ...) "DOU without X0 and with positive diagonal selection strength matrix (H1)."
#' @export
PCMParentClasses.DOU__noX0__posdiagH1 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noX0__posdiagH1 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec[!sapply(spec, is.null)]
}


#' @export
PCMDescribe.DOU__noX0__posdiagH2 <- function(model, ...) "DOU without X0 and with positive diagonal decorrelation rate matrix (H2)."
#' @export
PCMParentClasses.DOU__noX0__posdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noX0__posdiagH2 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__noX0__posdiagH2__posdiagSigma_x <- function(model, ...) "DOU without X0 and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__noX0__posdiagH2__posdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noX0__posdiagH2__posdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__noX0__posdiagH2__posdiagSigmae_x <- function(model, ...) "DOU without X0 and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigmae_x (i.e. no non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__noX0__posdiagH2__posdiagSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noX0__posdiagH2__posdiagSigmae_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec$Sigmae_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                        type = c("matrix", "diag", "positive.diag"),
                        description = "Positive diagonal Choleski factor of the non-phylogenetic variance-covariance matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__noX0__posdiagH2__posdiagSigma_x__posdiagSigmae_x <- function(model, ...) "DOU without X0 and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x and Sigmae_x (i.e. no phylogenetic or non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__noX0__posdiagH2__posdiagSigma_x__posdiagSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noX0__posdiagH2__posdiagSigma_x__posdiagSigmae_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec$Sigmae_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                        type = c("matrix", "diag", "positive.diag"),
                        description = "Positive diagonal Choleski factor of the non-phylogenetic variance-covariance matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__noX0__zeroH1__posdiagH2 <- function(model, ...) "DOU without X0 and with zero selection strength matrix (H1) and positive diagonal decorrelation rate matrix (H2)."
#' @export
PCMParentClasses.DOU__noX0__zeroH1__posdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noX0__zeroH1__posdiagH2 <- function(model, ...) {
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
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__noX0__zeroH1__posdiagH2__posdiagSigma_x <- function(model, ...) "DOU without X0 and with zero selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__noX0__zeroH1__posdiagH2__posdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noX0__zeroH1__posdiagH2__posdiagSigma_x <- function(model, ...) {
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
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__noX0__zeroH1__posdiagH2__posdiagSigmae_x <- function(model, ...) "DOU without X0 and with zero selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigmae_x (i.e. no non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__noX0__zeroH1__posdiagH2__posdiagSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noX0__zeroH1__posdiagH2__posdiagSigmae_x <- function(model, ...) {
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
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec$Sigmae_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                        type = c("matrix", "diag", "positive.diag"),
                        description = "Positive diagonal Choleski factor of the non-phylogenetic variance-covariance matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__noX0__zeroH1__posdiagH2__posdiagSigma_x__posdiagSigmae_x <- function(model, ...) "DOU without X0 and with zero selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x and Sigmae_x (i.e. no phylogenetic or non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__noX0__zeroH1__posdiagH2__posdiagSigma_x__posdiagSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noX0__zeroH1__posdiagH2__posdiagSigma_x__posdiagSigmae_x <- function(model, ...) {
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
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec$Sigmae_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                        type = c("matrix", "diag", "positive.diag"),
                        description = "Positive diagonal Choleski factor of the non-phylogenetic variance-covariance matrix")
  spec[!sapply(spec, is.null)]
}


#' @export
PCMDescribe.DOU__noX0__posdiagH1__posdiagH2 <- function(model, ...) "DOU without X0 and with positive diagonal selection strength matrix (H1) and positive diagonal decorrelation rate matrix (H2)."
#' @export
PCMParentClasses.DOU__noX0__posdiagH1__posdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noX0__posdiagH1__posdiagH2 <- function(model, ...) {
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
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__noX0__posdiagH1__posdiagH2__posdiagSigma_x <- function(model, ...) "DOU without X0 and with positive diagonal selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__noX0__posdiagH1__posdiagH2__posdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noX0__posdiagH1__posdiagH2__posdiagSigma_x <- function(model, ...) {
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
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__noX0__posdiagH1__posdiagH2__posdiagSigmae_x <- function(model, ...) "DOU without X0 and with positive diagonal selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigmae_x (i.e. no non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__noX0__posdiagH1__posdiagH2__posdiagSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noX0__posdiagH1__posdiagH2__posdiagSigmae_x <- function(model, ...) {
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
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec$Sigmae_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                        type = c("matrix", "diag", "positive.diag"),
                        description = "Positive diagonal Choleski factor of the non-phylogenetic variance-covariance matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__noX0__posdiagH1__posdiagH2__posdiagSigma_x__posdiagSigmae_x <- function(model, ...) "DOU without X0 and with positive diagonal selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x and Sigmae_x (i.e. no phylogenetic or non-phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__noX0__posdiagH1__posdiagH2__posdiagSigma_x__posdiagSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noX0__posdiagH1__posdiagH2__posdiagSigma_x__posdiagSigmae_x <- function(model, ...) {
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
                  description = "Non-negative diagonal decorrelation rate matrix")
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
PCMDescribe.DOU__noSigmae_x <- function(model, ...) "DOU without Sigmae_x."
#' @export
PCMParentClasses.DOU__noSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noSigmae_x <- function(model, ...) {
  spec <- NextMethod()
  spec$Sigmae_x <- NULL
  spec[!sapply(spec, is.null)]
}


#' @export
PCMDescribe.DOU__noSigmae_x__zeroH1 <- function(model, ...) "DOU without Sigmae_x and with zero selection strength matrix (H1)."
#' @export
PCMParentClasses.DOU__noSigmae_x__zeroH1 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noSigmae_x__zeroH1 <- function(model, ...) {
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
PCMDescribe.DOU__noSigmae_x__posdiagH1 <- function(model, ...) "DOU without Sigmae_x and with positive diagonal selection strength matrix (H1)."
#' @export
PCMParentClasses.DOU__noSigmae_x__posdiagH1 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noSigmae_x__posdiagH1 <- function(model, ...) {
  spec <- NextMethod()
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec[!sapply(spec, is.null)]
}


#' @export
PCMDescribe.DOU__noSigmae_x__posdiagH2 <- function(model, ...) "DOU without Sigmae_x and with positive diagonal decorrelation rate matrix (H2)."
#' @export
PCMParentClasses.DOU__noSigmae_x__posdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noSigmae_x__posdiagH2 <- function(model, ...) {
  spec <- NextMethod()
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__noSigmae_x__posdiagH2__posdiagSigma_x <- function(model, ...) "DOU without Sigmae_x and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__noSigmae_x__posdiagH2__posdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noSigmae_x__posdiagH2__posdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec[!sapply(spec, is.null)]
}


#' @export
PCMDescribe.DOU__noSigmae_x__zeroH1__posdiagH2 <- function(model, ...) "DOU without Sigmae_x and with zero selection strength matrix (H1) and positive diagonal decorrelation rate matrix (H2)."
#' @export
PCMParentClasses.DOU__noSigmae_x__zeroH1__posdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noSigmae_x__zeroH1__posdiagH2 <- function(model, ...) {
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
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__noSigmae_x__zeroH1__posdiagH2__posdiagSigma_x <- function(model, ...) "DOU without Sigmae_x and with zero selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__noSigmae_x__zeroH1__posdiagH2__posdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noSigmae_x__zeroH1__posdiagH2__posdiagSigma_x <- function(model, ...) {
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
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__noSigmae_x__posdiagH1__posdiagH2 <- function(model, ...) "DOU without Sigmae_x and with positive diagonal selection strength matrix (H1) and positive diagonal decorrelation rate matrix (H2)."
#' @export
PCMParentClasses.DOU__noSigmae_x__posdiagH1__posdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noSigmae_x__posdiagH1__posdiagH2 <- function(model, ...) {
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
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__noSigmae_x__posdiagH1__posdiagH2__posdiagSigma_x <- function(model, ...) "DOU without Sigmae_x and with positive diagonal selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__noSigmae_x__posdiagH1__posdiagH2__posdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noSigmae_x__posdiagH1__posdiagH2__posdiagSigma_x <- function(model, ...) {
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
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec[!sapply(spec, is.null)]
}




#######################
### noX0__noSigmae_x

#' @export
PCMDescribe.DOU__noX0__noSigmae_x <- function(model, ...) "DOU without X0 and Sigmae_x."
#' @export
PCMParentClasses.DOU__noX0__noSigmae_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noX0__noSigmae_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL
  spec[!sapply(spec, is.null)]
}


#' @export
PCMDescribe.DOU__noX0__noSigmae_x__zeroH1 <- function(model, ...) "DOU without X0 and without Sigmae_x and with zero selection strength matrix (H1)."
#' @export
PCMParentClasses.DOU__noX0__noSigmae_x__zeroH1 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noX0__noSigmae_x__zeroH1 <- function(model, ...) {
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
PCMDescribe.DOU__noX0__noSigmae_x__posdiagH1 <- function(model, ...) "DOU without X0 and without Sigmae_x and with positive diagonal selection strength matrix (H1)."
#' @export
PCMParentClasses.DOU__noX0__noSigmae_x__posdiagH1 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noX0__noSigmae_x__posdiagH1 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H1 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec[!sapply(spec, is.null)]
}


#' @export
PCMDescribe.DOU__noX0__noSigmae_x__posdiagH2 <- function(model, ...) "DOU without X0 and without Sigmae_x and with positive diagonal decorrelation rate matrix (H2)."
#' @export
PCMParentClasses.DOU__noX0__noSigmae_x__posdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noX0__noSigmae_x__posdiagH2 <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__noX0__noSigmae_x__posdiagH2__posdiagSigma_x <- function(model, ...) "DOU without X0 and without Sigmae_x and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__noX0__noSigmae_x__posdiagH2__posdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noX0__noSigmae_x__posdiagH2__posdiagSigma_x <- function(model, ...) {
  spec <- NextMethod()
  spec$X0 <- NULL
  spec$Sigmae_x <- NULL

  k <- attr(model, "k")
  regimes <- attr(model, "regimes")
  R <- length(regimes)

  spec$H2 <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                  type = c("matrix", "diag", "positive.diag"),
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec[!sapply(spec, is.null)]
}


#' @export
PCMDescribe.DOU__noX0__noSigmae_x__zeroH1__posdiagH2 <- function(model, ...) "DOU without X0 and without Sigmae_x and with zero selection strength matrix (H1) and positive diagonal decorrelation rate matrix (H2)."
#' @export
PCMParentClasses.DOU__noX0__noSigmae_x__zeroH1__posdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noX0__noSigmae_x__zeroH1__posdiagH2 <- function(model, ...) {
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
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__noX0__noSigmae_x__zeroH1__posdiagH2__posdiagSigma_x <- function(model, ...) "DOU without X0 and without Sigmae_x and with zero selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__noX0__noSigmae_x__zeroH1__posdiagH2__posdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noX0__noSigmae_x__zeroH1__posdiagH2__posdiagSigma_x <- function(model, ...) {
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
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__noX0__noSigmae_x__posdiagH1__posdiagH2 <- function(model, ...) "DOU without X0 and without Sigmae_x and with positive diagonal selection strength matrix (H1) and positive diagonal decorrelation rate matrix (H2)."
#' @export
PCMParentClasses.DOU__noX0__noSigmae_x__posdiagH1__posdiagH2 <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noX0__noSigmae_x__posdiagH1__posdiagH2 <- function(model, ...) {
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
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec[!sapply(spec, is.null)]
}

#' @export
PCMDescribe.DOU__noX0__noSigmae_x__posdiagH1__posdiagH2__posdiagSigma_x <- function(model, ...) "DOU without X0 and without Sigmae_x and with positive diagonal selection strength matrix (H1) and with positive diagonal decorrelation rate matrix (H2) and positive diagonal Sigma_x (i.e. no phylogenetic correlation between the traits)."
#' @export
PCMParentClasses.DOU__noX0__noSigmae_x__posdiagH1__posdiagH2__posdiagSigma_x <- function(model) c("DOU", "GaussianPCM", "PCM")
#' @export
PCMSpecifyParams.DOU__noX0__noSigmae_x__posdiagH1__posdiagH2__posdiagSigma_x <- function(model, ...) {
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
                  description = "Non-negative diagonal decorrelation rate matrix")
  spec$Sigma_x <- list(default = array(0, dim = c(k, k, R), dimnames = list(NULL, NULL, regimes)),
                       type = c("matrix", "diag", "positive.diag"),
                       description = "Positive diagonal Choleski factor of the unit-time variance-covariance matrix of the BM-process")
  spec[!sapply(spec, is.null)]
}
