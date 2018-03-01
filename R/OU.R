#' Validate OU parameters,
#' @export
PCMValidate.OU <- function(tree, model, verbose = FALSE) {
  if(verbose) {
    print('Validating model...')
  }

  if(is.null(model$Sigmae) |
     is.null(dim(model$Sigmae)) |
     length(dim(model$Sigmae)) != 3) {
    stop("ERR:02201:PCMBase:OU.R:PCMValidate.OU:: Expecting the model to have a member called Sigmae with dimensions k x k x R, where R is the number of regimes and k is the number of traits.")
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
      tree = tree, modelName = "OU",
      k = k, R = R, regimesUnique = regimesUnique,
      paramNames = list("H", "Theta", "Sigma", "Sigmae"),
      paramStorageModes = list("double", "double", "double", "double"),
      paramDims = list(c(k, k, R), c(k, R), c(k, k, R), c(k, k, R))
    ),
    verbose = verbose)
}

#' Create a conditional multivariate OU distribution
#'
#' @param H,Theta,Sigma parameters of the multivariate OU process; H is a k x k
#' matrix, Theta is a k-vector and Sigma is a k x k matrix
#'
#' @return a list containging the passed parameters as well as
#' a function `random` of arguments n (number of observation k-vectors to generate),
#' x0 (initial k-vector of values), t (numeric time); and a function `density` for
#' calculating the density of multivariate vector under the specified distribution
#' and given an initial value and time.
#'
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @importFrom expm expm
#'
#' @export
PCMCond.OU <- function(tree, model, r=1, verbose=FALSE) {
  with(model, {
    H <- as.matrix(model$H[,,r])
    Theta <- model$Theta[,r]
    Sigma <- as.matrix(model$Sigma[,,r])

    if(length(unique(c(length(Theta), dim(H), dim(Sigma)))) != 1) {
      stop('ERR:02221:PCMBase:OU.R:PCMCond.OU:: Some of H, Theta or Sigma has a wrong dimension.')
    }

    V <- PCMCondVOU(H, Sigma)
    omega <- function(t, edgeIndex, e_Ht = NULL) {
      if(is.null(e_Ht)) {
        e_Ht <- expm(-t*H)
      }
      I <- diag(nrow(H))
      (I-e_Ht) %*% Theta
    }
    Phi <- function(t, edgeIndex, e_Ht = NULL) {
      if(is.null(e_Ht)) {
        expm(-t*H)
      } else {
        e_Ht
      }
    }
    random <- function(n=1, x0, t, edgeIndex) {
      e_Ht <- expm(-t*H)
      rmvnorm(n=n, mean = omega(t, edgeIndex, e_Ht) + Phi(t, edgeIndex, e_Ht) %*% x0, sigma=V(t, edgeIndex))
    }
    density <- function(x, x0, t, edgeIndex, log=FALSE) {
      e_Ht <- expm(-t*H)
      dmvnorm(x, mean = omega(t, edgeIndex, e_Ht) + Phi(t, edgeIndex, e_Ht)%*%x0, sigma=V(t, edgeIndex), log=log)
    }

    list(omega = omega, Phi = Phi, V = V, random=random, density=density)
  })
}

