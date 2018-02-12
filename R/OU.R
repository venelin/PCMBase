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

PLambdaP_1.OU <- function(H) {
  # here the argument H is a matrix specifying the alphas in a OU process
  r <- eigen(H)
  if(isTRUE(all.equal(rcond(r$vectors),0))) {
    stop(paste0("ERR:02211:PCMBase:OU.R:PLambdaP_1.OU:: The provided H matrix is defective. Its matrix of eigenvectors is computationally singular. H=", toString(H)))
  }
  list(lambda=r$values, P=r$vectors, P_1=solve(r$vectors))
}

#' sums of pairs of elements in a vector
#' @param lambda
Lambda_ij.OU <- function(lambda) {
  sapply(lambda, function(li) sapply(lambda, function(lj) li+lj))
}

#' Create a function of time that calculates (1-exp(-lambda_{ij}*time))/lambda_{ij}
#' for every element lambda_{ij} of the input matrix Lambda_ij
#' @param Lambda_ij a squared numerical matrix of dimension k x k
#' @param threshold.Lambda_ij a 0-threshold for abs(Lambda_i + Lambda_j), where Lambda_i and Lambda_j are
#'   eigenvalues of the parameter matrix H. This threshold-values is used as a condition to
#'   take the limit time of the expression `(1-exp(-Lambda_ij*time))/Lambda_ij` as
#'   `(Lambda_i+Lambda_j) --> 0`. You can control this value by the global option
#'   "PCMBase.Threshold.Lambda_ij". The default value (1e-8) is suitable for branch lengths bigger than
#'   1e-6. For smaller branch lengths, you are may want to increase the threshold value using, e.g.
#'   `options(PCMBase.Threshold.Lambda_ij=1e-6)`.
#' @return a function of time returning a matrix with entries formed from the
#'  above function or the limit time if |Lambda_{ij}|<=trehshold0.
fLambda_ij.OU <- function(
  Lambda_ij,
  threshold.Lambda_ij = getOption("PCMBase.Threshold.Lambda_ij", 1e-8)) {

  idx0 <- which(abs(Lambda_ij)<=threshold.Lambda_ij)
  function(time) {
    res <- (1-exp(-Lambda_ij*time))/Lambda_ij
    if(length(idx0)>0) {
      res[idx0] <- time
    }
    res
  }
}


#' Generate a multivariate (MV) OU variance-covariance function
#' @param lambda vector of the eignevalues of the matrix H of a MV OU process.
#' @param P matrix of the eigenvectors of the matrix H
#' @param P_1 inverse eigenvectors matrix
#' @param Sigma the matrix Sigma of a MV OU process.
#' @param Sigmae optional additional error (environmental) variance covariance
#' matrix or single number (if the error has the same variance for all dimensions).
#' This matrix is added to the MV-OU variance covariance matrix.
#' @param threshold.Lambda_ij a 0-threshold for abs(Lambda_i + Lambda_j), where Lambda_i and Lambda_j are
#'   eigenvalues of the parameter matrix H. This threshold-values is used as a condition to
#'   take the limit time of the expression `(1-exp(-Lambda_ij*time))/Lambda_ij` as
#'   `(Lambda_i+Lambda_j) --> 0`. You can control this value by the global option
#'   "PCMBase.Threshold.Lambda_ij". The default value (1e-8) is suitable for branch lengths bigger than
#'   1e-6. For smaller branch lengths, you are may want to increase the threshold value using, e.g.
#'   `options(PCMBase.Threshold.Lambda_ij=1e-6)`.
#' @return a function of one numerical argument (time), which calculates the
#' expected variance covariance matrix of a MV-OU process after time, given
#' the specified arguments.
V.OU <- function(
  lambda, P, P_1, Sigma,
  threshold.Lambda_ij = getOption("PCMBase.Threshold.Lambda_ij", 1e-8) ) {

  Lambda_ij <- Lambda_ij.OU(lambda)
  fLambda_ij <- fLambda_ij.OU(Lambda_ij, threshold.Lambda_ij)

  P_1SigmaP_t <- P_1 %*% Sigma %*% t(P_1)

  # need to evoque P as well to make it available for the daughter function
  force(P)
  force(P_1)

  function(time) {
    Re(P %*% (fLambda_ij(time)*P_1SigmaP_t) %*% t(P))
  }
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
  PLP_1 <- PLambdaP_1.OU(H)
  fV <- V.OU(PLP_1$lambda, PLP_1$P, PLP_1$P_1, Sigma)

  random <- function(n=1, x0, t, e) {
    e_Ht <- expm(-t*H)
    I <- diag(nrow(H))
    rmvnorm(n=n, mean=e_Ht%*%x0 + (I-e_Ht) %*% Theta, sigma=fV(t))
  }
  density <- function(x, x0, t, e, log=FALSE) {
    e_Ht <- expm(-t*H)
    I <- diag(nrow(H))
    dmvnorm(x, mean=e_Ht%*%x0 + (I-e_Ht) %*% Theta, sigma=fV(t), log=log)
  }

  list(H=H, Theta=Theta, Sigma=Sigma, random=random, density=density, vcov=fV)
  })
}

#' Calculate the coefficients A, b, C, d, E, f of the general
#' form (eq. 1) for each edge in a tree
#'
#' @param tree a phylo object (see details)
#' @param model parameters of the OU process. This must be a
#' named list with the following elements:
#' H: a k x k x R array, where R is the number of regimes of the
#' OU process, k is the number of variables (traits), each H[,,r]
#' containing the matrix H for regime r;
#' Theta: a k x R matrix, row Theta[, r] containing the long-term
#' mean Theta for regime r;
#' Sigma: a k x k x R array, each Sigma[,,r] containing the
#' matrix Sigma for regime r;
#' Sigmae: a k x k x R array, each Sigmae[,,r] representing a diagonal matrix
#' with elements on the diagona corresponding to the environmental variances for
#' the k traits in regime r
#' @param pc a k x M logical matrix representing the present coordinates at each
#' node
#'
#' @details The dimnames
#'
#' @return a named list containing the following elements:
#' A: a k x k x M array, A[,,i] corresponding to Ai for
#' each branch ending at node or tip i;
#' b: a k x M matrix, b[,i] corresponding to the vectors bi;
#' C: a k x k x M array, C[,,i] corresponding to the
#' matrices Ci;
#' d: a k x M matrix, d[,i] corresponding to the vectors di;
#' E: a k x k x M array, E[,,i] corresponding to the matrices Ei;
#' f: a vector, f[i] correspondign to fi
#'
#' @importFrom expm expm
#' @export
PCMAbCdEf.OU <- function(tree, model,
                      metaI=PCMValidate.OU(tree, model, verbose=verbose),
                      pc, verbose=FALSE) {
  # number of regimes
  R <- metaI$R

  # number of tips
  N <- metaI$N

  # number of traits (variables)
  k <- metaI$k

  # number of nodes
  M <- metaI$M

  tree <- tree

  P <- array(NA, dim=c(k, k, R), dimnames=dimnames(model$H))
  P_1 <- array(NA, dim=c(k, k, R), dimnames=dimnames(model$H))
  lambda <- array(NA, dim=c(k, R), dimnames=dimnames(model$H)[-2])

  fV.OU <- list()

  for(r in 1:R) {
    PLambdaP_1 <- PLambdaP_1.OU(model$H[,,r])
    P[,,r] <- PLambdaP_1$P
    P_1[,,r] <- PLambdaP_1$P_1
    lambda[,r] <- PLambdaP_1$lambda

    # create the V.OU function for regime r
    fV.OU[[r]] <- V.OU(lambda[,r], as.matrix(P[,,r]), as.matrix(P_1[,,r]),
                       as.matrix(model$Sigma[,,r]))
  }

  V <- array(NA, dim=c(k, k, M))
  V_1 <- array(NA, dim=c(k, k, M))
  e_Ht <- array(NA, dim=c(k, k, M))

  # returned general form parameters
  A <- array(NA, dim=c(k, k, M))
  b <- array(NA, dim=c(k, M))
  C <- array(NA, dim=c(k, k, M))
  d <- array(NA, dim=c(k, M))
  E <- array(NA, dim=c(k, k, M))
  f <- array(NA, dim=c(M))

  # vector of regime indices for each branch
  r <- metaI$regimes

  # identity k x k matrix
  I <- diag(k)

  # iterate over the edges
  for(e in 1:(M-1)) {
    # parent node
    j <- tree$edge[e, 1]
    # daughter node
    i <- tree$edge[e, 2]

    # length of edge leading to i
    ti <- tree$edge.length[e]

    # present coordinates in parent and daughter nodes
    kj <- pc[,j]
    ki <- pc[,i]

    V[,,i] <- fV.OU[[r[e]]](ti)

    if(i <= N) {
      # add environmental variance at each tip node
      V[,,i] <- V[,,i] + model$Sigmae[,,r[e]]
    }

    V_1[ki,ki,i] <- solve(V[ki,ki,i])
    e_Ht[,,i] <- expm(-ti*as.matrix(model$H[,,r[e]]))

    # now compute PCMAbCdEf according to eq (8) in doc.
    # here A is from the general form (not the alpa from OU process)
    A[ki,ki,i] <- -0.5*V_1[ki,ki,i]

    b[ki,i] <- V_1[ki,ki,i] %*% (I[ki,]-e_Ht[ki,,i]) %*% model$Theta[,r[e]]

    C[kj,kj,i] <- -0.5*t(matrix(e_Ht[ki,kj,i], sum(ki), sum(kj))) %*% V_1[ki,ki,i] %*% e_Ht[ki,kj,i]

    d[kj,i] <- -t(matrix(e_Ht[ki,kj,i], sum(ki), sum(kj))) %*% V_1[ki,ki,i] %*% (matrix(I[ki,]-e_Ht[ki,,i], sum(ki))) %*% model$Theta[,r[e]]

    E[kj,ki,i] <- t(matrix(e_Ht[ki,kj,i], sum(ki), sum(kj))) %*% V_1[ki,ki,i]

    f[i] <-
      -0.5*(sum(ki)*log(2*pi) + log(det(as.matrix(V[ki,ki,i]))) +
              t(model$Theta[,r[e]]) %*% t(matrix(I[ki,]-e_Ht[ki,,i], sum(ki))) %*%
              V_1[ki,ki,i] %*% (matrix(I[ki,]-e_Ht[ki,,i], sum(ki))) %*% model$Theta[,r[e]])
  }

  list(A=A, b=b, C=C, d=d, E=E, f=f, V=V, V_1=V_1)
}
