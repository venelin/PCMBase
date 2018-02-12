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


#' Generate a multivariate (MV) BM variance-covariance function
#' @param Sigma the matrix Sigma of a MV BM process.
#' @param threshold0
#' @return a function of one numerical argument (time), which calculates the
#' expected variance covariance matrix of a MV-BM process after time, given
#' the specified arguments.

V.BM <- function(Sigma, threshold0=0) {

  force(Sigma)
  fun = function(time) {
        return (time*Sigma)
  }
  return (fun)
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

    fV <- V.BM(Sigma)

    random <- function(n=1, x0, t, e) {
      rmvnorm(n=n,
              mean=x0,
              sigma=fV(t))
    }
    density <- function(x, x0, t, e, log=FALSE) {
      dmvnorm(x,
              mean=x0,
              sigma=fV(t), log=log)
    }

    list(Sigma=Sigma, random=random, density=density, vcov=fV)
  })
}



#' Calculate the coefficients A, b, C, d, E, f of the general
#' form (eq. 1) for each edge in a tree
#'
#' @param tree a phylo object (see details)
#' @param model parameters of the BM process. This must be a
#' named list with the following elements:
#' Sigma: a k x k x R array, each Sigma[,,r] containing the
#' matrix Sigma for regime r;
#' Sigmae: a k x k x R array, each Sigmae[,,r] representing a diagonal matrix
#' with elements on the diagona corresponding to the environmental variances for
#' the k traits in regime r
#' @param presentCoords a k x M logical matrix representing the present coordinates at each
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
#' @export
PCMAbCdEf.BM <- function(tree, model,
                      metaI=PCMValidate.BM(tree, model, verbose=verbose),
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
    fV.BM[[r]] <- V.BM(as.matrix(model$Sigma[,,r]))

  }

  V <- array(NA, dim=c(k, k, M))
  V_1 <- array(NA, dim=c(k, k, M))

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

    # present coordinates in parent and daughte nodes
    kj <- pc[,j]
    ki <- pc[,i]

    V[,,i] <- fV.BM[[r[e]]](ti)


    if(i<=N) {
      # add environmental variance at each tip node
      V[,,i] <- V[,,i] + model$Sigmae[,,r[e]]
    }

    V_1[ki,ki,i] <- solve(V[ki,ki,i])

    # now compute PCMAbCdEf according to eq (16) in doc.
    # here A is from the general form
    A[ki,ki,i] <- (-0.5*V_1[ki,ki,i])

    b[ki,i] <- 0

    C[kj,kj,i] <- (-0.5*(t(matrix(I[ki,kj], sum(ki), sum(kj))) %*% V_1[ki,ki,i]) %*% matrix(I[ki,kj], sum(ki), sum(kj)))

    d[kj,i] <- 0

    E[kj,ki,i] <- (t(matrix(I[ki,kj], sum(ki), sum(kj)))%*%V_1[ki,ki,i])

    f[i] <- -0.5*(sum(ki)*log(2*pi) + log(det(as.matrix(V[ki,ki,i]))))
  }

  list(A=A, b=b, C=C, d=d, E=E, f=f, V=V, V_1=V_1)
}
