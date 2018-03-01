#' Validate BM parameters
#' @export
PCMValidate.BM <- function(tree, model, verbose=FALSE) {
  if(verbose) {
    print('Validating model...')
  }
  if(is.null(model$Sigmae) | is.null(dim(model$Sigmae))) {
    stop("ERR:02101:PCMBase:BM.R:PCMValidate.BM:: Expecting the model to have a member called Sigmae with dimensions k x k x R, where R is the number of regimes and k is the number of traits.")
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
      tree = tree, modelName = "BM",
      k = k, R = R, regimesUnique = regimesUnique,
      paramNames = list("Sigma", "Sigmae"),
      paramStorageModes = list("double", "double"),
      paramDims = list(c(k, k, R), c(k, k, R))
    ),
    verbose = verbose)
}

#' Create a conditional multivariate BM distribution
#' @param Sigma of the multivariate BM process; Sigma is a k x k matrix
#' @return a list containging the passed parameters as well as
#' a function `random` of arguments n (number of observation k-vectors to generate),
#' x0 (initial k-vector of values), t (numeric time); and a function `density` for
#' calculating the density of multivariate vector under the specified distribution
#' and given an initial value and time.
#' @importFrom mvtnorm rmvnorm dmvnorm
#'
#' @export
PCMCond.BM <- function(tree, model, r=1, verbose=FALSE) {
  with(model, {
    Sigma <- as.matrix(model$Sigma[,,r])

    if(length(unique(c(dim(Sigma))))!=1) {
      # this is a dummy check to evaluate Theta
      print(paste('dim(Sigma)=', dim(Sigma)))
      stop('ERR:02102:PCMBase:BM.R:PCMCond.BM:: Sigma has a wrong dimension.')
    }

    V <- PCMCondVOU(matrix(0, nrow(Sigma), ncol(Sigma)), Sigma)
    omega <- function(t, edgeIndex) {
      rep(0, nrow(Sigma))
    }
    Phi <- function(t, edgeIndex, e_Ht = NULL) {
      diag(nrow(Sigma))
    }
    random <- function(n=1, x0, t, edgeIndex) {
      rmvnorm(n=n, mean=x0, sigma=V(t))
    }
    density <- function(x, x0, t, edgeIndex, log=FALSE) {
      dmvnorm(x, mean=x0, sigma=V(t), log=log)
    }

    list(omega = omega, Phi = Phi, V = V, random=random, density=density)
  })
}
