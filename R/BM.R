#' Validate BM parameters
#'
validateModel.BM <- function(tree, model, verbose=FALSE) {
  if(verbose) {
    print('Validating model...')
  }
  if(!all(c('Sigma', 'Sigmae')%in%names(model))) {
    stop("model should be a named list with elements 'Sigma', 'Sigmae'")
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

  # number of traits
  k <- dim(model$Sigma)[2]
  if(verbose) {
    cat('k=', k, '\n')
  }

  if(!is.array(model$Sigma)|!is.array(model$Sigmae) |
     !isTRUE(all.equal(dim(model$Sigma[regimes,,,drop=FALSE]), c(R, k, k))) |
     !isTRUE(all.equal(dim(model$Sigmae[regimes,,,drop=FALSE]), c(R, k, k)))) {
    if(verbose) {
      print('dim Sigma:')
      print(dim(model$Sigma))
      print('dim Sigmae:')
      print(dim(model$Sigmae))
    }
    stop("Incorrect dimensions for some of the parameters Sigma and Sigmae")
  }

  if(R>1 & (!all.equal(dimnames(model$Sigma)[[1]], dimnames(model$Sigmae)[[1]]))) {
    if(verbose) {
      print('dimnames Sigma:')
      print(dimnames(model$Sigma))
      print('dimnames Sigmae:')
      print(dimnames(model$Sigmae))
    }
    stop("Disagreeing regime-names (dimnames(X)[[1]]) of Sigma and Sigmae")
  }

  if(R==1) {
    regimes <- rep(1, length(tree$edge.length))
  } else {
    regimes <- match(tree$edge.regime, dimnames(model$Sigma)[[1]])
    if(any(is.na(regimes))) {
      if(verbose) {
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


#' Generate a multivariate (MV) BM variance-covariance function
#' @param Sigma the matrix Sigma of a MV BM process.
#' @param threshold0
#' @return Sigma.

V.BM <- function(Sigma, threshold0=0) {

  dummy = Sigma*2

  fun = function(time) {
        return (time*Sigma)
  }
  return (fun)
}

#' Create a conditional multivariate BM distribution
#' @param Sigma of the multivariate BM process; Sigma is a k x k matrix
#' @return a list containging the passed parameters as well as
#' a function mvr of arguments n (number of observation k-vectors to generate),
#' x0 (initial k-vector of values), t (numeric time); and a function mvd for
#' calculating the density of multivariate vector under the specified distribution
#' and given an initial value and time.
mvcond.BM <- function(model, r=1, verbose=FALSE) {
  with(model, {
    Sigma <- as.matrix(model$Sigma[r,,])

    if(length(unique(c(dim(Sigma))))!=1) {
      # this is a dummy check to evaluate Theta
      print(paste('dim(Sigma)=', dim(Sigma)))
      stop('Sigma has a wrong dimension.')
    }

    fV <- V.BM(Sigma)

    mvr <- function(n=1, x0, t) {
      mvtnorm::rmvnorm(n=n,
                       mean=x0,
                       sigma=fV(t))
    }
    mvd <- function(x, x0, t, log=FALSE) {
      dmvnorm(x,
              mean=x0,
              sigma=fV(t), log=log)
    }

    list(Sigma=Sigma, mvr=mvr, mvd=mvd, vcov=fV)
  })
}



#' Calculate the coefficients A, b, C, d, E, f of the general
#' form (eq. 1) for each edge in a tree
#'
#' @param tree a phylo object (see details)
#' @param model parameters of the BM process. This must be a
#' named list with the following elements:
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
AbCdEf.BM <- function(tree, model,
                      metaI=validateModel.BM(tree, model, verbose=verbose),
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

  fV.BM <- list()

  for(r in 1:R) {

    # create the V.BM function for regime r
    fV.BM[[r]] <- V.BM(as.matrix(model$Sigma[r,,]))

  }

  V <- array(NA, dim=c(M, k, k))
  V_1 <- array(NA, dim=c(M, k, k))

  # returned general form parameters
  A <- array(NA, dim=c(M, k, k))
  b <- array(NA, dim=c(M, k))
  C <- array(NA, dim=c(M, k, k))
  d <- array(NA, dim=c(M, k))
  E <- array(NA, dim=c(M, k, k))
  f <- array(NA, dim=c(M))

  # vector of regime indices for each branch
  r <- metaI$regimes

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

    #print(fV.BM[[r[e]]](ti))

    V[i,,] <- fV.BM[[r[e]]](ti)


    if(i<=N) {
      # add environmental variance at each tip node
      V[i,,] <- V[i,,] + model$Sigmae[r[e],,]
    }

    V_1[i,ki,ki] <- solve(V[i,ki,ki])

    # now compute AbCdEf according to eq (16) in doc.
    # here A is from the general form
    A[i,ki,ki] <- (-0.5*V_1[i,ki,ki])

    b[i,ki] <- 0

    C[i,kj,kj] <- (-0.5*V_1[i,ki,ki])

    d[i,kj] <- 0

    E[i,kj,ki] <- (V_1[i,ki,ki])

    f[i] <- -0.5*(sum(ki)*log(2*pi) + log(det(as.matrix(V[i,ki,ki]))))
  }

  list(A=A, b=b, C=C, d=d, E=E, f=f,V=V)
}
