#' Validate TwoSpeedOU parameters
#' @export
PCMValidate.TwoSpeedOU <- function(tree, model, verbose = FALSE) {
  if(verbose) {
    cat('Validating model...')
  }
  if(is.null(model$Sigmae) | is.null(dim(model$Sigmae))) {
    stop("ERR:02401:PCMBase:TwoSpeedOU.R:PCMValidate.TwoSpeedOU:: Expecting the model to have a member called Sigmae with dimensions k x k x R, where R is the number of regimes and k is the number of traits.")
  }

  R <- dim(model$Sigmae)[3]
  k <- dim(model$Sigmae)[1]
  regimesUnique <- dimnames(model$Sigmae)[[3]]

  if(is.null(regimesUnique)) {
    regimesUnique <- 1:dim(model$Sigmae)[[3]]
  }

  PCMValidateGeneral(
    tree = tree, model = model,
    modelSpec = PCMSpecify(
      tree = tree, modelName = "TwoSpeedOU",
      k = k, R = R, regimesUnique = regimesUnique,
      paramNames = list("H1", "H2", "Theta", "Sigma", "Sigmae"),
      paramStorageModes = list("double", "double", "double", "double", "double"),
      paramDims = list(c(k, k, R), c(k, k, R), c(k, R), c(k, k, R), c(k, k, R))
    ),
    verbose = verbose)
}

#' Create a conditional multivariate TwoSpeedOU distribution
#' @param H1,H2,Theta,Sigma parameters of the multivariate TwoSpeedOU process; H1,H2 are
#'  k x k matrices, Theta is a k-vector and Sigma is a k x k matrix
#' @return a list containging the passed parameters as well as
#' a function `random` of arguments n (number of observation k-vectors to generate),
#' x0 (initial k-vector of values), t (numeric time); and a function `density` for
#' calculating the density of multivariate vector under the specified distribution
#' and given an initial value and time.
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @importFrom expm expm
#'
#' @export
PCMCond.TwoSpeedOU <- function(tree, model, r=1, verbose=FALSE) {
  with(model, {
    H1 <- as.matrix(model$H1[,,r])
    H2 <- as.matrix(model$H2[,,r])
    Theta <- model$Theta[,r]
    Sigma <- as.matrix(model$Sigma[,,r])

    if(length(unique(c(length(Theta), dim(H1), dim(H2), dim(Sigma))))!=1) {
      stop('ERR:02421:PCMBase:TwoSpeedOU.R:PCMCond.TwoSpeedOU:: Some of H1, H2, Theta or Sigma has a wrong dimension.')
    }

    V <- PCMCondVOU(H2, Sigma)
    omega <- function(t, edgeIndex, e_H1t = NULL) {
      if(is.null(e_H1t)) {
        e_H1t <- expm(-t*H1)
      }
      I <- diag(nrow(H1))
      (I-e_H1t) %*% Theta
    }
    Phi <- function(t, edgeIndex, e_H1t = NULL) {
      if(is.null(e_H1t)) {
        expm(-t*H1)
      } else {
        e_H1t
      }
    }
    random <- function(n=1, x0, t, edgeIndex) {
      e_H1t <- expm(-t*H1)
      rmvnorm(n=n, mean = omega(t, edgeIndex, e_H1t) + Phi(t, edgeIndex, e_H1t)%*%x0, sigma=V(t, edgeIndex))
    }
    density <- function(x, x0, t, edgeIndex, log=FALSE) {
      e_H1t <- expm(-t*H1)
      dmvnorm(x, mean=omega(t, edgeIndex, e_H1t) + Phi(t, edgeIndex, e_H1t)%*%x0, sigma=V(t, edgeIndex), log=log)
    }
    list(omega = omega, Phi = Phi, V = V, random=random, density=density)
  })
}
