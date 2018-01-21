#' Validate JOU parameters
#' @export
validateModel.JOU <- function(tree, model, verbose=FALSE) {
  if(verbose) {
    print('Validating model...')
  }
  if(is.null(model$Sigmae) | is.null(dim(model$Sigmae))) {
    stop("ERR:02301:PCMBase:JOU.R:validateModel.JOU:: Expecting the model to have a member called Sigmae with dimensions
         R x k x k, where R is the number of regimes and k is the number of
         traits.")
  }

  R <- dim(model$Sigmae)[1]
  k <- dim(model$Sigmae)[2]
  regimesUnique <- dimnames(model$Sigmae)[[1]]

  if(is.null(regimesUnique)) {
    regimesUnique <- 1:dim(model$Sigmae)[[1]]
  }

  if(is.null(tree$edge.jump)) {
    stop("ERR:02302:PCMBase:JOU.R:validateModel.JOU:: Expecting the tree to have a member edge.jump - an integer vector of
         0's and 1's describing if there is a jump for each branch of the tree.")
  }

  if(!all(tree$edge.jump %in% as.integer(0:1))) {
    stop("ERR:02303:PCMBase:JOU.R:validateModel.JOU:: Check that tree$edge.jump is an integer vector of 0's and 1's")
  }

  if(length(tree$edge.jump) != nrow(tree$edge)) {
    stop("ERR:02304:PCMBase:JOU.R:validateModel.JOU:: Check that tree$edge.jump has nrow(tree$edge) elements.")
  }

  validateModelGeneral(
    tree = tree, model = model,
    modelSpec = specifyModel(
      tree = tree, modelName = "JOU",
      k = k, R = R, regimesUnique = regimesUnique,
      paramNames = list("Alpha", "Theta", "Sigma", "Sigmae", "mj", "Sigmaj"),
      paramStorageModes = list("double", "double", "double", "double", "double", "double"),
      paramDims = list(c(R, k, k), c(R, k), c(R, k, k), c(R, k, k), c(R, k), c(R, k, k))
    ),
    verbose = verbose)
}

PLambdaP_1.JOU <- function(Alpha) {
  # here the argument Alpha is an Alpha matrix specifying the alphas in a JOU process
  r <- eigen(Alpha)
  if(isTRUE(all.equal(rcond(r$vectors),0))) {
    stop(paste0("ERR:02311:PCMBase:JOU.R:PLambdaP_1.JOU:: The provided Alpha matrix is defective. Its matrix of eigenvectors is computationally singular. Alpha=", toString(Alpha)))
  }
  list(lambda=r$values, P=r$vectors, P_1=solve(r$vectors))
}

#' sums of pairs of elements in a vector
#' @param lambda
Lambda_ij.JOU <- function(lambda) {
  sapply(lambda, function(li) sapply(lambda, function(lj) li+lj))
}

#' Create a function of time that calculates (1-exp(-lambda_{ij}*time))/lambda_{ij}
#' for every element lambad_{ij} of the input matrix Lambda_ij
#' @param Lambda_ij a squared numerical matrix of dimension k x k
#' @param threshold.Lambda_ij a 0-threshold for abs(Lambda_i + Lambda_j), where Lambda_i and Lambda_j are
#'   eigenvalues of the parameter matrix Alpha. This threshold-values is used as a condition to
#'   take the limit time of the expression `(1-exp(-Lambda_ij*time))/Lambda_ij` as
#'   `(Lambda_i+Lambda_j) --> 0`. You can control this value by the global option
#'   "PCMBase.Threshold.Lambda_ij". The default value of 1e-8 is suitable for branch lengths bigger than
#'   1e-6. For smaller branch lengths, you are likely to increase the threshold value using, e.g.
#'   `options(PCMBase.Threshold.Lambda_ij = 1e-6)`.
#' @return a function of time returning a matrix with entries formed from the
#'  above function or the limit time if |Lambda_{ij}|<=trehshold0.
fLambda_ij.JOU <- function(
  Lambda_ij,
  threshold.Lambda_ij = getOption("PCMBase.Threshold.Lambda_ij", 1e-8)) {

  idx0 <- which(abs(Lambda_ij) <= threshold.Lambda_ij)
  function(time) {
    res <- (1-exp(-Lambda_ij*time))/Lambda_ij
    if(length(idx0)>0) {
      res[idx0] <- time
    }
    res
  }
}

#' Create a function of time that calculates exp(-Alpha*time)
#' @param the argument Alpha is an Alpha matrix kxk specifying the alphas in a JOU process
#'
#' @importFrom expm expm
#' @return a function of time returning a matrix
texp.JOU <- function(Alpha){
  force(Alpha)
  return(function(time){
    res = expm(-time*Alpha)
    return(res)
  })
}

#' Generate a multivariate (MV) JOU variance-covariance function
#' @param lambda vector of the eignevalues of the matrix Alpha of a MV JOU process.
#' @param P matrix of the eigenvectors of the matrix Alpha
#' @param P_1 inverse eigenvectors matrix
#' @param Sigma the matrix Sigma of a MV JOU process.
#' @param Alpha is the matrix Alpha of a MV JOU process.
#' @param Sigmaj is the variance matrix of the jump distribution in JOU process
#' @param threshold.Lambda_ij a 0-threshold for abs(Lambda_i + Lambda_j), where Lambda_i and Lambda_j are
#'   eigenvalues of the parameter matrix Alpha. This threshold-values is used as a condition to
#'   take the limit time of the expression `(1-exp(-Lambda_ij*time))/Lambda_ij` as
#'   `(Lambda_i+Lambda_j) --> 0`. You can control this value by the global option
#'   "PCMBase.Threshold.Lambda_ij". The default value (1e-8) is suitable for branch lengths bigger than
#'   1e-6. For smaller branch lengths, you are may want to increase the threshold value using, e.g.
#'   `options(PCMBase.Threshold.Lambda_ij=1e-6)`.
#' @return a function of two numerical arguments (time) and (xi:binary value denoting jump or not) , which calculates the
#' expected variance covariance matrix of a MV-JOU process after time, given
#' the specified arguments.
#' @export
V.JOU <- function(
  lambda, P, P_1, Sigma, Alpha, Sigmaj,
  threshold.Lambda_ij = getOption("PCMBase.Threshold.Lambda_ij", 1e-8)) {

  Lambda_ij <- Lambda_ij.JOU(lambda)
  fLambda_ij <- fLambda_ij.JOU(Lambda_ij, threshold.Lambda_ij)
  P_1SigmaP_t <- P_1%*%Sigma%*%t(P_1)
  force(Alpha)
  force(Sigmaj)
  e_At <- texp.JOU(Alpha)

  force(P)
  force(P_1)

  function(time, xi) {
    e_Ati <- e_At(time)
    Re((P %*% (fLambda_ij(time)*P_1SigmaP_t) %*% t(P)) + xi*(e_Ati %*% Sigmaj %*% t(e_Ati)))
  }
}

#' Create a conditional multivariate JOU distribution
#' @param Alpha,Theta,Sigma,Sigmaj,mj,xi parameters of the multivariate JOU process; Alpha is a k x k
#' matrix, Theta,mj are k-vectors and Sigma,Sigmaj are k x k matrices,xi is a vector of length equal to the number
#' of edges in the tree
#' @return a list containging the passed parameters as well as
#' a function mvr of arguments n (number of observation k-vectors to generate),
#' x0 (initial k-vector of values), t (numeric time); and a function mvd for
#' calculating the density of multivariate vector under the specified distribution
#' and given an initial value and time.
#' @importFrom expm expm
#' @importFrom mvtnorm rmvnorm dmvnorm
mvcond.JOU <- function(tree, model, r=1, verbose=FALSE) {
  with(model, {
    Alpha <- as.matrix(model$Alpha[r,,])
    Theta <- model$Theta[r,]
    Sigma <- as.matrix(model$Sigma[r,,])
    Sigmaj <- as.matrix(model$Sigmaj[r,,])
    mj <- model$mj[r,]
    xi <- tree$edge.jump

    if(length(unique(c(length(Theta), dim(Alpha), dim(mj), dim(Sigmaj), dim(Sigma))))!=1) {
      stop('ERR:02321:PCMBase:JOU.R:mvcond.JOU:: Some of Alpha, Theta, Sigma,  Sigmaj or mj have a wrong dimension.')
    }
    PLP_1 <- PLambdaP_1.JOU(Alpha)
    fV <- V.JOU(PLP_1$lambda, PLP_1$P, PLP_1$P_1, Sigma, Alpha, Sigmaj)

    mvr <- function(n=1, x0, t, e) {
      e_At <- expm(-t*Alpha)
      I <- diag(nrow(Alpha))
      rmvnorm(n=n,
              mean=e_At%*%x0 + (I-e_At)%*%Theta + mj*xi[e],
              sigma=fV(t, xi[e]))
    }
    mvd <- function(x, x0, t, e, log=FALSE) {
      e_At <- expm(-t*Alpha)
      I <- diag(nrow(Alpha))
      dmvnorm(x,
              mean=e_At%*%x0 + (I-e_At)%*%Theta + mj*xi[e],
              sigma=fV(t, xi[e]), log=log)
    }

    list(Alpha=Alpha, Theta=Theta, Sigma=Sigma, Sigmaj = Sigmaj, mj = mj, mvr=mvr, mvd=mvd, vcov=fV)
  })
}



#' Calculate the coefficients A, b, C, d, E, f of the general
#' form (eq. 1) for each edge in a tree
#'
#' @param tree a phylo object (see details)
#' @param model parameters of the JOU process. This must be a
#' named list with the following elements:
#' Alpha: a R x k x k array, where R is the number of regimes of the
#' JOU process, k is the number of variables (traits), each Alpha[r,,]
#' containing the matrix Alpha for regime r;
#' Theta: a R x k matrix, row Theta[r, ] containing the long-term
#' mean Theta for regime r;
#' Sigma: a R x k x k array, each Sigma[r,,] containing the
#' matrix Sigma for regime r;
#' Sigmae: a R x k x k array, each Sigmae[r,,] representing a diagonal matrix
#' with elements on the diagona corresponding to the environmental variances for
#' the k traits in regime r
#' Sigmaj: a R x k x k array, each Sigmaj[r,,] containing the
#' matrix Sigmaj for regime r;
#' mj: a R x k matrix, row mj[r, ] containing mj for regime r;
#' xi: a M-1 binary vector indicating jump or not
#' @param pc a M x k logical matrix representing the present coordinates at each
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
#' @importFrom expm expm
#'
#' @export
AbCdEf.JOU <- function(tree, model,
                       metaI=validateModel.JOU(tree, model, verbose=verbose),
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
  P <- array(NA, dim=c(R, k, k), dimnames=dimnames(model$Alpha))
  P_1 <- array(NA, dim=c(R, k, k), dimnames=dimnames(model$Alpha))
  lambda <- array(NA, dim=c(R, k), dimnames=dimnames(model$Alpha)[-3])

  fV.JOU <- list()

  for(r in 1:R) {
    PLambdaP_1 <- PLambdaP_1.JOU(model$Alpha[r,,])
    P[r,,] <- PLambdaP_1$P
    P_1[r,,] <- PLambdaP_1$P_1
    lambda[r,] <- PLambdaP_1$lambda

    fV.JOU[[r]] <- V.JOU(lambda[r,],
                         as.matrix(P[r,,]),
                         as.matrix(P_1[r,,]),
                         as.matrix(model$Sigma[r,,]),
                         as.matrix(model$Alpha[r,,]),
                         as.matrix(model$Sigmaj[r,,]))
  }

  V <- array(NA, dim=c(M, k, k))
  V_1 <- array(NA, dim=c(M, k, k))
  e_At <- array(NA, dim=c(M, k, k))
  e_ATt <- array(NA, dim=c(M, k, k))

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
    # binary vector indicating jump or not per edge
    xi <- tree$edge.jump
    # present coordinates in parent and daughte nodes
    kj <- pc[j,]
    ki <- pc[i,]
    V[i,,] <- fV.JOU[[r[e]]](ti, xi[e])

    if(i<=N) {
      # add environmental variance at each tip node
      V[i,,] <- V[i,,] + model$Sigmae[r[e],,]
    }

    V_1[i,ki,ki] <- solve(V[i,ki,ki])
    e_At[i,,] <- expm(-ti*as.matrix(model$A[r[e],,]))

    # now compute AbCdEf
    # here A is from the general form (not the alpha from JOU process)
    A[i,ki,ki] <- -0.5*V_1[i,ki,ki]

    b[i,ki] <- V_1[i,ki,ki] %*%
      ((matrix(I[ki,]-e_At[i,ki,], sum(ki))) %*% model$Theta[r[e],] + e_At[i,ki,] %*% model$mj[r[e],]*xi[e])

    C[i,kj,kj] <- -0.5*t(matrix(e_At[i,ki,kj], sum(ki), sum(kj))) %*% V_1[i,ki,ki] %*% matrix(e_At[i,ki,kj], sum(ki), sum(kj))

    d[i,kj] <- -t(matrix(e_At[i,ki,kj], sum(ki), sum(kj))) %*% V_1[i,ki,ki] %*%
      ((matrix(I[ki,]-e_At[i,ki,], sum(ki))) %*% model$Theta[r[e],] + e_At[i,ki,] %*% model$mj[r[e],] * xi[e])

    E[i,kj,ki] <- t(matrix(e_At[i,ki,kj], sum(ki), sum(kj))) %*% V_1[i,ki,ki]

    f[i] <-
      -0.5*(sum(ki)*log(2*pi) + log(det(as.matrix(V[i,ki,ki]))) +
              (t(model$Theta[r[e],]) %*% t(matrix(I[ki,]-e_At[i,ki,], sum(ki))) + t(model$mj[r[e],]) %*% t(matrix(e_At[i,ki,], sum(ki)))*xi[e]) %*%
              V_1[i,ki,ki] %*% ((matrix(I[ki,]-e_At[i,ki,], sum(ki))) %*% model$Theta[r[e],] + e_At[i,ki,] %*% model$mj[r[e],]*xi[e]))
  }

  list(A=A, b=b, C=C, d=d, E=E, f=f, e_At=e_At, V=V, V_1=V_1)
}

