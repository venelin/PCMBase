#' Validate OU parameters
#'
validateModel.OU <- function(tree, model, verbose=FALSE) {
  if(verbose) {
    print('Validating model...')
  }
  if(!all(c('Alpha', 'Theta', 'Sigma', 'Sigmae')%in%names(model))) {
    stop("model should be a named list with elements 'Alpha', 'Theta', 'Sigma', 'Sigmae'")
  }

  # number of tips
  N <- length(tree$tip.label)

  # number of nodes including the root
  M <- nrow(tree$edge)+1
  if(verbose) {
    cat('M=', M, '\n')
  }

  # number of regimes
  if(is.null(tree$edge.regime)) {
    regimes <- 1
    R <- 1
  } else {
    regimes <- unique(tree$edge.regime)
    R <- length(regimes)
  }

  if(verbose) {
    cat('R=',R,'\n')
  }

  if(!is.matrix(model$Theta) | nrow(model$Theta)<R) {
    stop(paste0('Theta should be a R x k matrix, where R is the number of regimes and k - the number of traits: R=',R))
  }

  # number of traits
  k <- dim(model$Theta)[2]
  if(verbose) {
    cat('k=', k, '\n')
  }

  if(!is.array(model$Alpha)|!is.array(model$Sigma)|!is.array(model$Sigmae) |
     !all.equal(dim(model$Alpha[regimes,,,drop=FALSE]),c(R, k, k)) |
     !all.equal(dim(model$Sigma[regimes,,,drop=FALSE]), c(R, k, k)) |
     !all.equal(dim(model$Sigmae[regimes,,,drop=FALSE]), c(R, k, k))) {
    if(verbose) {
      print('dim Alpha:')
      print(dim(model$Alpha))
      print('dim Sigma:')
      print(dim(model$Sigma))
      print('dim Sigmae:')
      print(dim(model$Sigmae))
    }
    stop("Incorrect dimensions for some of the parameters Alpha, Sigma and Sigmae")
  }

  if(R>1 & (!all.equal(dimnames(model$Alpha)[[1]], dimnames(model$Theta)[[1]]) |
            !all.equal(dimnames(model$Alpha)[[1]], dimnames(model$Sigma)[[1]]) |
            !all.equal(dimnames(model$Alpha)[[1]], dimnames(model$Sigmae)[[1]]))) {
    if(verbose) {
      print('dimnames Alpha:')
      print(dimnames(model$Alpha))
      print('dimnames Theta:')
      print(dimnames(model$Theta))
      print('dimnames Sigma:')
      print(dimnames(model$Sigma))
      print('dimnames Sigmae:')
      print(dimnames(model$Sigmae))
    }
    stop("Disagreeing regime-names (dimnames(X)[[1]]) of Alpha, Theta, Sigma and Sigmae")
  }

  if(R==1) {
    regimes <- rep(1, length(tree$edge.length))
  } else {
    regimes <- match(tree$edge.regime, dimnames(model$Alpha)[[1]])
    if(any(is.na(regimes))) {
      if(verbose) {
        print('dimnames Alpha:')
        print(dimnames(model$Alpha))
        print('dimnames Theta:')
        print(dimnames(model$Theta))
        print('dimnames Sigma:')
        print(dimnames(model$Sigma))
        print('dimnames Sigmae:')
        print(dimnames(model$Sigmae))
      }
      stop("Not all regime-names (names in tree$edge.length) are specified parameters")
    }
  }
  list(N=N, M=M, R=R, k=k, regimes=regimes)
}

PLambdaP_1.OU <- function(Alpha) {
  # here the argument Alpha is an Alpha matrix specifying the alphas in a OU process
  r <- eigen(Alpha)
  list(lambda=r$values, P=r$vectors, P_1=solve(r$vectors))
}

#' sums of pairs of elements in a vector
#' @param lambda
Lambda_ij.OU <- function(lambda) {
  sapply(lambda, function(li) sapply(lambda, function(lj) li+lj))
}

#' Create a function of time that calculates (1-exp(-lambda_{ij}*time))/lambda_{ij}
#' for every element lambad_{ij} of the input matrix Lambda_ij
#' @param Lambda_ij a squared numerical matrix of dimension k x k
#' @param threshold0 numeric threshold value. For |Lambda_{ij}| smaller or equal
#'  than threshold the limit time is returned.
#' @return a function of time returning a matrix with entries formed from the
#'  above function or the limit time if |Lambda_{ij}|<=trehshold0.
fLambda_ij.OU <- function(Lambda_ij, threshold0=0) {
  idx0 <- which(abs(Lambda_ij)<=threshold0)
  function(time) {
    res <- (1-exp(-Lambda_ij*time))/Lambda_ij
    if(length(idx0)>0) {
      res[idx0] <- time
    }
    res
  }
}


#' Generate a multivariate (MV) OU variance-covariance function
#' @param lambda vector of the eignevalues of the matrix Alpha of a MV OU process.
#' @param P matrix of the eigenvectors of the matrix Alpha
#' @param P_1 inverse eigenvectors matrix
#' @param Sigma the matrix Sigma of a MV OU process.
#' @param Sigmae optional additional error (environmental) variance covariance
#' matrix or single number (if the error has the same variance for all dimensions).
#' This matrix is added to the MV-OU variance covariance matrix.
#' @param threshold0
#' @return a function of one numerical argument (time), which calculates the
#' expected variance covariance matrix of a MV-OU process after time, given
#' the specified arguments.
V.OU <- function(lambda, P, P_1, Sigma, threshold0=0) {
  Lambda_ij <- Lambda_ij.OU(lambda)
  fLambda_ij <- fLambda_ij.OU(Lambda_ij, threshold0)
  P_1SigmaP_t <- P_1%*%Sigma%*%t(P_1)

  # need to evoque P as well to make it available for the daughter function
  if(!all.equal(dim(P),dim(P_1))) {
    # Dummy code: this should never happen
    stop('Error!')
  }

  function(time) {
    P %*% (fLambda_ij(time)*P_1SigmaP_t) %*% t(P)
  }
}

#' Create a conditional multivariate OU distribution
#' @param Alpha,Theta,Sigma parameters of the multivariate OU process; Alpha is a k x k
#' matrix, Theta is a k-vector and Sigma is a k x k matrix
#' @return a list containging the passed parameters as well as
#' a function mvr of arguments n (number of observation k-vectors to generate),
#' x0 (initial k-vector of values), t (numeric time); and a function mvd for
#' calculating the density of multivariate vector under the specified distribution
#' and given an initial value and time.
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @importFrom expm expm
#'
mvcond.OU <- function(model, r=1, verbose=FALSE) {
  with(model, {
    Alpha <- as.matrix(model$Alpha[r,,])
    Theta <- model$Theta[r,]
    Sigma <- as.matrix(model$Sigma[r,,])

  if(length(unique(c(length(Theta), dim(Alpha), dim(Sigma))))!=1) {
    # this is a dummy check to evaluate Theta
    print(paste('dim(Alpha)=', dim(Alpha)))
    print(paste('length(Theta)=', length(Theta)))
    print(paste('dim(Sigma)=', dim(Sigma)))
    stop('Some of Alpha, Theta or Sigma has a wrong dimension.')
  }
  PLP_1 <- PLambdaP_1.OU(Alpha)
  fV <- V.OU(PLP_1$lambda, PLP_1$P, PLP_1$P_1, Sigma)

  mvr <- function(n=1, x0, t, e) {
    e_At <- expm::expm(-t*Alpha)
    I <- diag(nrow(Alpha))
    rmvnorm(n=n, mean=e_At%*%x0 + (I-e_At) %*% Theta, sigma=fV(t))
  }
  mvd <- function(x, x0, t, e, log=FALSE) {
    e_At <- expm::expm(-t*Alpha)
    I <- diag(nrow(Alpha))
    dmvnorm(x, mean=e_At%*%x0 + (I-e_At) %*% Theta, sigma=fV(t), log=log)
  }

  list(Alpha=Alpha, Theta=Theta, Sigma=Sigma, mvr=mvr, mvd=mvd, vcov=fV)
  })
}



#' Calculate the coefficients A, b, C, d, E, f of the general
#' form (eq. 1) for each edge in a tree
#'
#' @param tree a phylo object (see details)
#' @param model parameters of the OU process. This must be a
#' named list with the following elements:
#' Alpha: a R x k x k array, where R is the number of regimes of the
#' OU process, k is the number of variables (traits), each Alpha[r,,]
#' containing the matrix Alpha for regime r;
#' Theta: a R x k matrix, row Theta[r, ] containing the long-term
#' mean Theta for regime r;
#' Sigma: a R x k x k array, each Sigma[r,,] containing the
#' matrix Sigma for regime r;
#' Sigmae: a R x k x k array, each Sigmae[r,,] representing a diagonal matrix
#' with elements on the diagona corresponding to the environmental variances for
#' the k traits in regime r
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
AbCdEf.OU <- function(tree, model,
                      metaI=validateModel.OU(tree, model, verbose=verbose),
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

  fV.OU <- list()

  for(r in 1:R) {
    PLambdaP_1 <- PLambdaP_1.OU(model$Alpha[r,,])
    P[r,,] <- PLambdaP_1$P
    P_1[r,,] <- PLambdaP_1$P_1
    lambda[r,] <- PLambdaP_1$lambda

    # create the V.OU function for regime r
    fV.OU[[r]] <- V.OU(lambda[r,], as.matrix(P[r,,]), as.matrix(P_1[r,,]),
                       as.matrix(model$Sigma[r,,]))
  }

  V <- array(NA, dim=c(M, k, k))
  V_1 <- array(NA, dim=c(M, k, k))
  e_At <- array(NA, dim=c(M, k, k))

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

    V[i,,] <- fV.OU[[r[e]]](ti)

    if(i<=N) {
      # add environmental variance at each tip node
      V[i,,] <- V[i,,] + model$Sigmae[r[e],,]
    }

    V_1[i,ki,ki] <- solve(V[i,ki,ki])
    e_At[i,,] <- expm::expm(-ti*as.matrix(model$A[r[e],,]))

    # now compute AbCdEf according to eq (8) in doc.
    # here A is from the general form (not the alpa from OU process)
    A[i,ki,ki] <- -0.5*V_1[i,ki,ki]

    b[i,ki] <- V_1[i,ki,ki] %*% (I[ki,]-e_At[i,ki,]) %*% model$Theta[r[e],]

    C[i,kj,kj] <- -0.5*t(e_At[i,ki,kj]) %*% V_1[i,ki,ki] %*% e_At[i,ki,kj]

    d[i,kj] <- -t(e_At[i,ki,kj]) %*% V_1[i,ki,ki] %*% (I[ki,]-e_At[i,ki,]) %*% model$Theta[r[e],]

    E[i,kj,ki] <- t(e_At[i,ki,kj]) %*% V_1[i,ki,ki]

    f[i] <-
      -0.5*(sum(ki)*log(2*pi) + log(det(as.matrix(V[i,ki,ki]))) +
              t(model$Theta[r[e],]) %*% t(I[,ki]-e_At[i,,ki]) %*%
              V_1[i,ki,ki] %*% (I[,ki]-e_At[i,,ki]) %*% model$Theta[r[e],])
  }

  list(A=A, b=b, C=C, d=d, E=E, f=f, e_At=e_At, V=V)
}
