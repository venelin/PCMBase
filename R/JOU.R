#' Validate JOU parameters
#' @export
PCMValidate.JOU <- function(tree, model, verbose=FALSE) {
  if(verbose) {
    print('Validating model...')
  }
  if(is.null(model$Sigmae) | is.null(dim(model$Sigmae))) {
    stop("ERR:02301:PCMBase:JOU.R:PCMValidate.JOU:: Expecting the model to have a member called Sigmae with dimensions k x k x R, where R is the number of regimes and k is the number of traits.")
  }

  R <- dim(model$Sigmae)[3]
  k <- dim(model$Sigmae)[1]
  regimesUnique <- dimnames(model$Sigmae)[[3]]

  if(is.null(regimesUnique)) {
    regimesUnique <- 1:dim(model$Sigmae)[[3]]
  }

  if(is.null(tree$edge.jump)) {
    stop("ERR:02302:PCMBase:JOU.R:PCMValidate.JOU:: Expecting the tree to have a member edge.jump - an integer vector of 0's and 1's describing if there is a jump for each branch of the tree.")
  }

  if(!all(tree$edge.jump %in% as.integer(0:1))) {
    stop("ERR:02303:PCMBase:JOU.R:PCMValidate.JOU:: Check that tree$edge.jump is an integer vector of 0's and 1's")
  }

  if(length(tree$edge.jump) != nrow(tree$edge)) {
    stop("ERR:02304:PCMBase:JOU.R:PCMValidate.JOU:: Check that tree$edge.jump has nrow(tree$edge) elements.")
  }

  PCMValidateGeneral(
    tree = tree, model = model,
    modelSpec = PCMSpecify(
      tree = tree, modelName = "JOU",
      k = k, R = R, regimesUnique = regimesUnique,
      paramNames = list("H", "Theta", "Sigma", "Sigmae", "mj", "Sigmaj"),
      paramStorageModes = list("double", "double", "double", "double", "double", "double"),
      paramDims = list(c(k, k, R), c(k, R), c(k, k, R), c(k, k, R), c(k, R), c(k, k, R))
    ),
    verbose = verbose)
}

#' Create a conditional multivariate JOU distribution
#' @param H,Theta,Sigma,Sigmaj,mj,xi parameters of the multivariate JOU process; H is a k x k
#' matrix, Theta,mj are k-vectors and Sigma,Sigmaj are k x k matrices,xi is a vector of length equal to the number
#' of edges in the tree
#' @return a list containging the passed parameters as well as
#' a function `random` of arguments n (number of observation k-vectors to generate),
#' x0 (initial k-vector of values), t (numeric time); and a function `density` for
#' calculating the density of multivariate vector under the specified distribution
#' and given an initial value and time.
#' @importFrom expm expm
#' @importFrom mvtnorm rmvnorm dmvnorm
#'
#' @export
PCMCond.JOU <- function(tree, model, r=1, verbose=FALSE) {
  with(model, {
    H <- as.matrix(model$H[,,r])
    Theta <- model$Theta[,r]
    Sigma <- as.matrix(model$Sigma[,,r])
    Sigmaj <- as.matrix(model$Sigmaj[,,r])
    mj <- model$mj[,r]
    xi <- tree$edge.jump

    if(length(unique(c(length(Theta), dim(H), dim(mj), dim(Sigmaj), dim(Sigma))))!=1) {
      stop('ERR:02321:PCMBase:JOU.R:PCMCond.JOU:: Some of H, Theta, Sigma,  Sigmaj or mj have a wrong dimension.')
    }

    V <- PCMCondVOU(H, Sigma, Sigmaj, xi)
    omega <- function(t, edgeIndex, e_Ht = NULL) {
      if(is.null(e_Ht)) {
        e_Ht <- expm(-t*H)
      }
      I <- diag(nrow(H))
      xi[edgeIndex] * e_Ht%*%mj +  (I-e_Ht)%*%Theta
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
      rmvnorm(n=n, mean = omega(t, edgeIndex, e_Ht) + Phi(t, edgeIndex, e_Ht) %*% x0, sigma=V(t, edgeIndex, e_Ht))
    }
    density <- function(x, x0, t, edgeIndex, log=FALSE) {
      e_Ht <- expm(-t*H)
      dmvnorm(x, mean = omega(t, edgeIndex, e_Ht) + Phi(t, edgeIndex, e_Ht) %*% x0, sigma=V(t, edgeIndex, e_Ht), log=log)
    }

    list(omega = omega, Phi = Phi, V = V, random=random, density=density)
  })
}
