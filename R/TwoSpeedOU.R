#' Validate TwoSpeedOU parameters
#' @export
validateModel.TwoSpeedOU <- function(tree, model, verbose = FALSE) {
  if(verbose) {
    cat('Validating model...')
  }
  if(is.null(model$Sigmae) | is.null(dim(model$Sigmae))) {
    stop("ERR:02401:PCMBase:TwoSpeedOU.R:validateModel.TwoSpeedOU:: Expecting the model to have a member called Sigmae with dimensions R x k x k, where R is the number of regimes and k is the number of traits.")
  }

  R <- dim(model$Sigmae)[1]
  k <- dim(model$Sigmae)[2]
  regimesUnique <- dimnames(model$Sigmae)[[1]]

  if(is.null(regimesUnique)) {
    regimesUnique <- 1:dim(model$Sigmae)[[1]]
  }

  validateModelGeneral(
    tree = tree, model = model,
    modelSpec = specifyModel(
      tree = tree, modelName = "TwoSpeedOU",
      k = k, R = R, regimesUnique = regimesUnique,
      paramNames = list("Alpha1", "Alpha2", "Theta", "Sigma", "Sigmae"),
      paramStorageModes = list("double", "double", "double", "double", "double"),
      paramDims = list(c(R, k, k), c(R, k, k), c(R, k), c(R, k, k), c(R, k, k))
    ),
    verbose = verbose)
}

PLambdaP_1.TwoSpeedOU <- function(Alpha) {
  # here the argument Alpha is an Alpha matrix specifying the alphas in a TwoSpeedOU process
  r <- eigen(Alpha)
  if(isTRUE(all.equal(rcond(r$vectors),0))) {
    stop(paste0("ERR:02411:PCMBase:TwoSpeedOU.R:PLambdaP_1.TwoSpeedOU:: The provided Alpha matrix is defective - its matrix of eigenvectors is computationally singular. Alpha=", toString(Alpha)))
  }
  list(lambda=r$values, P=r$vectors, P_1=solve(r$vectors))
}

#' sums of pairs of elements in a vector
#' @param lambda
Lambda_ij.TwoSpeedOU <- function(lambda) {
  sapply(lambda, function(li) sapply(lambda, function(lj) li+lj))
}

#' Create a function of time that calculates (1-exp(-lambda_{ij}*time))/lambda_{ij}
#' for every element lambad_{ij} of the input matrix Lambda_ij
#' @param Lambda_ij a squared numerical matrix of dimension k x k
#' @param threshold.Lambda_ij a 0-threshold for abs(Lambda_i + Lambda_j), where Lambda_i
#'   and Lambda_j are eigenvalues of the parameter matrix Alpha. This
#'   threshold-value is used as a condition to
#'   take the limit time of the expression `(1-exp(-Lambda_ij*time))/Lambda_ij` as
#'   `(Lambda_i+Lambda_j) --> 0`. You can control this value by the global option
#'   "PCMBase.Threshold.Lambda_ij". The default value (1e-8) is suitable for branch
#'    lengths bigger than 1e-6. For smaller branch lengths, you are may want to
#'    increase the threshold value using, e.g.
#'   `options(PCMBase.Threshold.Lambda_ij=1e-6)`.
#' @return a function of time returning a matrix with entries formed from the
#'  above function or the limit time if |Lambda_{ij}|<=trehshold0.
fLambda_ij.TwoSpeedOU <- function(
  Lambda_ij,
  threshold.Lambda_ij = getOption("PCMBase.Threshold.Lambda_ij", 1e-8) ) {

  idx0 <- which(abs(Lambda_ij)<=threshold.Lambda_ij)
  function(time) {
    res <- (1-exp(-Lambda_ij*time))/Lambda_ij
    if(length(idx0)>0) {
      res[idx0] <- time
    }
    res
  }
}


#' Generate a multivariate (MV) TwoSpeedOU variance-covariance function
#' @param lambda vector of the eignevalues of the matrix Alpha2 of a MV TwoSpeedOU process.
#' @param P matrix of the eigenvectors of the matrix Alpha
#' @param P_1 inverse eigenvectors matrix
#' @param Sigma the matrix Sigma of a MV TwoSpeedOU process.
#' @param threshold.Lambda_ij a 0-threshold for abs(Lambda_i + Lambda_j), where Lambda_i
#'   and Lambda_j are eigenvalues of the parameter matrix Alpha. This
#'   threshold-value is used as a condition to
#'   take the limit time of the expression `(1-exp(-Lambda_ij*time))/Lambda_ij` as
#'   `(Lambda_i+Lambda_j) --> 0`. You can control this value by the global option
#'   "PCMBase.Threshold.Lambda_ij". The default value (1e-8) is suitable for branch
#'    lengths bigger than 1e-6. For smaller branch lengths, you are may want to
#'    increase the threshold value using, e.g.
#'   `options(PCMBase.Threshold.Lambda_ij=1e-6)`.
#' @return a function of one numerical argument (time), which calculates the
#' expected variance covariance matrix of a MV-TwoSpeedOU process after time, given
#' the specified arguments.
V.TwoSpeedOU <- function(
  lambda, P, P_1, Sigma,
  threshold.Lambda_ij = getOption("PCMBase.Threshold.Lambda_ij", 1e-8) ) {

  Lambda_ij <- Lambda_ij.TwoSpeedOU(lambda)
  fLambda_ij <- fLambda_ij.TwoSpeedOU(Lambda_ij, threshold.Lambda_ij)
  P_1SigmaP_t <- P_1%*%Sigma%*%t(P_1)

  force(P)
  force(P_1)

  function(time) {
    Re(P %*% (fLambda_ij(time)*P_1SigmaP_t) %*% t(P))
  }
}

#' Create a conditional multivariate TwoSpeedOU distribution
#' @param Alpha1,Alpha2,Theta,Sigma parameters of the multivariate TwoSpeedOU process; Alpha1,Alpha2 are
#'  k x k matrices, Theta is a k-vector and Sigma is a k x k matrix
#' @return a list containging the passed parameters as well as
#' a function mvr of arguments n (number of observation k-vectors to generate),
#' x0 (initial k-vector of values), t (numeric time); and a function mvd for
#' calculating the density of multivariate vector under the specified distribution
#' and given an initial value and time.
mvcond.TwoSpeedOU <- function(tree, model, r=1, verbose=FALSE) {
  with(model, {
    Alpha1 <- as.matrix(model$Alpha1[r,,])
    Alpha2 <- as.matrix(model$Alpha2[r,,])
    Theta <- model$Theta[r,]
    Sigma <- as.matrix(model$Sigma[r,,])

    if(length(unique(c(length(Theta), dim(Alpha1), dim(Alpha2), dim(Sigma))))!=1) {
      stop('ERR:02421:PCMBase:TwoSpeedOU.R:mvcond.TwoSpeedOU:: Some of Alpha1, Alpha2, Theta or Sigma has a wrong dimension.')
    }
    PLP_1 <- PLambdaP_1.TwoSpeedOU(Alpha2)
    fV <- V.TwoSpeedOU(PLP_1$lambda, PLP_1$P, PLP_1$P_1, Sigma)

    mvr <- function(n=1, x0, t, e) {
      e_A1t <- expm::expm(-t*Alpha1)
      I <- diag(nrow(Alpha1))
      mvtnorm::rmvnorm(n=n,
                       mean=e_A1t%*%x0 + (I-e_A1t)%*%Theta,
                       sigma=fV(t))
    }
    mvd <- function(x, x0, t, e, log=FALSE) {
      e_A1t <- expm::expm(-t*Alpha1)
      I <- diag(nrow(Alpha1))
      dmvnorm(x,
              mean=e_A1t%*%x0 + (I-e_A1t)%*%Theta,
              sigma=fV(t), log=log)
    }

    list(Alpha1=Alpha1, Alpha2 = Alpha2, Theta=Theta, Sigma=Sigma, mvr=mvr, mvd=mvd, vcov=fV)
  })
}



#' Calculate the coefficients A, b, C, d, E, f of the general
#' form (eq. 1) for each edge in a tree
#'
#' @param tree a phylo object (see details)
#' @param model parameters of the TwoSpeedOU process. This must be a
#' named list with the following elements:
#' Alpha1, Alpha2: R x k x k arrays, where R is the number of regimes of the
#' TwoSpeedOU process, k is the number of variables (traits)
#' Theta: a R x k matrix, row Theta[r, ] containing the long-term
#' mean Theta for regime r;
#' Sigma: a R x k x k array, each Sigma[r,,] containing the
#' matrix Sigma for regime r;
#' Sigmae: a R x k x k array, each Sigmae[r,,] representing a diagonal matrix
#' with elements on the diagona corresponding to the environmental variances for
#' the k traits in regime r
#' @param presentCoords a M x k logical matrix representing the present coordinates at each
#' node
#'
#' @details The dimnames
#'
#' @return a named list containing the following elements:
#' A: a M x k x k array, A[i,,] corresponding to Ai for
#' each branch ending at node or tip i;
#' b: a M x k matrix, b[i,] corresponding to the vectors bi;
#' C: a M x k x k array, C[i,,] corresponding to the
#' matrices Ci;
#' d: a M x k matrix, d[i,] corresponding to the vectors di;
#' E: a M x k x k array, E[i,,] corresponding to the matrices Ei;
#' f: a vector, f[i] correspondign to fi
#'
#' @export
AbCdEf.TwoSpeedOU <- function(tree, model,
                      metaI=validateModel.TwoSpeedOU(tree, model, verbose=verbose),
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

  P <- array(NA, dim=c(R, k, k), dimnames=dimnames(model$Alpha2))
  P_1 <- array(NA, dim=c(R, k, k), dimnames=dimnames(model$Alpha2))
  lambda <- array(NA, dim=c(R, k), dimnames=dimnames(model$Alpha2)[-3])

  fV.TwoSpeedOU <- list()

  for(r in 1:R) {
    PLambdaP_1 <- PLambdaP_1.TwoSpeedOU(model$Alpha2[r,,])
    P[r,,] <- PLambdaP_1$P
    P_1[r,,] <- PLambdaP_1$P_1
    lambda[r,] <- PLambdaP_1$lambda

    # create the V.TwoSpeedOU function for regime r
    fV.TwoSpeedOU[[r]] <- V.TwoSpeedOU(
      lambda[r,], as.matrix(P[r,,]), as.matrix(P_1[r,,]),
      as.matrix(model$Sigma[r,,]))
  }

  V <- array(NA, dim=c(M, k, k))
  V_1 <- array(NA, dim=c(M, k, k))
  e_A1t <- array(NA, dim=c(M, k, k))

  # returned general form parameters
  A <- array(NA, dim=c(M, k, k))
  b <- array(NA, dim=c(M, k))
  C <- array(NA, dim=c(M, k, k))
  d <- array(NA, dim=c(M, k))
  E <- array(NA, dim=c(M, k, k))
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

    # present coordinates in parent and daughte nodes
    kj <- pc[j,]
    ki <- pc[i,]

    V[i,,] <- fV.TwoSpeedOU[[r[e]]](ti)

    if(i<=N) {
      # add environmental variance at each tip node
      V[i,,] <- V[i,,] + model$Sigmae[r[e],,]
    }

    V_1[i,ki,ki] <- solve(V[i,ki,ki])
    e_A1t[i,,] <- expm::expm(-ti*as.matrix(model$Alpha1[r[e],,]))

    # now compute AbCdEf
    # here A is from the general form
    A[i,ki,ki] <- -0.5*V_1[i,ki,ki]

    b[i,ki] <- V_1[i,ki,ki] %*% (I[ki,]-e_A1t[i,ki,]) %*% model$Theta[r[e],]

    C[i,kj,kj] <- -0.5*t(matrix(e_A1t[i,ki,kj], sum(ki), sum(kj))) %*% V_1[i,ki,ki] %*% e_A1t[i,ki,kj]

    d[i,kj] <- -t(matrix(e_A1t[i,ki,kj], sum(ki), sum(kj))) %*% V_1[i,ki,ki] %*% (I[ki,]-e_A1t[i,ki,]) %*% model$Theta[r[e],]

    E[i,kj,ki] <- t(matrix(e_A1t[i,ki,kj], sum(ki), sum(kj))) %*% V_1[i,ki,ki]

    f[i] <-
      -0.5*(sum(ki)*log(2*pi) + log(det(as.matrix(V[i,ki,ki]))) +
              t(model$Theta[r[e],]) %*% t(matrix(I[ki,]-e_A1t[i,ki,], sum(ki))) %*%
              V_1[i,ki,ki] %*% (matrix(I[ki,]-e_A1t[i,ki,], sum(ki))) %*% model$Theta[r[e],])
  }

  list(A=A, b=b, C=C, d=d, E=E, f=f, e_At=e_A1t, V=V, V_1=V_1)
}
