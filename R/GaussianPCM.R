# Copyright 2018 Venelin Mitov
#
# This file is part of PCMBase.
#
# PCMBase is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PCMBase is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PCMBase.  If not, see <http://www.gnu.org/licenses/>.

#' @name PCMBase
#'
#' @title A General Framework for Gaussian Phylogenetic Comparative Models
#'
#' @description
NULL

#' @export
PCMParentClasses.GaussianPCM <- function(model) {
  "PCM"
}

#' @inherit PCMCond
#' @return For GaussianPCM models, a named list with the following members:
#' \item{omega}{d}
#' \item{Phi}{}
#' \item{V}{}
#' @export
PCMCond.GaussianPCM <- function(tree, model, r=1, metaI=PCMInfo(NULL, tree, model, verbose), verbose = FALSE) {
  stop(paste('ERR:02111:PCMBase:GaussianPCM.R:PCMCond.GaussianPCM:: This is an abstract function which only defines an interface but should not be called explicitly. Possibly you forgot implementing PCMCond for a daughter class.'))
}

#' @export
PCMSim.GaussianPCM <- function(
  tree, model, X0,
  metaI = PCMInfo(X = NULL, tree = tree, model = model, verbose = verbose),
  verbose = FALSE) {

  if(length(X0)!=metaI$k) {
    stop(paste('ERR:02102:PCMBase:GaussianPCM.R:PCMSim:: X0 must be of length', metaI$k, '.'))
  }

  values <- matrix(0, nrow=metaI$k, ncol=dim(tree$edge)[1]+1)
  values[, metaI$N + 1] <- X0

  ordBF <- metaI$preorder

  # create a list of random generator functions for each regime
  PCMCondObjects <- lapply(1:metaI$RModel, function(r) {
    PCMCond(tree, model = model, r = r, metaI = metaI, verbose = verbose)
  })

  for(edgeIndex in ordBF) {
    obj <- PCMCondObjects[[metaI$r[edgeIndex]]]
    if(!is.null(obj$random)) {
      values[, tree$edge[edgeIndex, 2]] <-
        obj$random(n=1, x0 = values[, tree$edge[edgeIndex,1]], t = tree$edge.length[edgeIndex], edgeIndex = edgeIndex)
    } else {
      values[, tree$edge[edgeIndex, 2]] <-
        PCMCondRandom(obj, n=1, x0 = values[, tree$edge[edgeIndex,1]], t = tree$edge.length[edgeIndex], edgeIndex = edgeIndex, metaI = metaI)
    }
  }

  values
}

#' @export
PCMLik.GaussianPCM <- function(
  X, tree, model,
  metaI = PCMInfo(X, tree, model, verbose = verbose),
  log = TRUE,
  verbose = FALSE) {

  # will change this value if there is no error
  value.NA <- getOption("PCMBase.Value.NA", as.double(NA))

  PCMLmr <- PCMLmr(X, tree, model, metaI, verbose = verbose, root.only = TRUE)#,
            # silent = TRUE)

  if(class(PCMLmr) == "try-error") {
    errL <- PCMParseErrorMessage(PCMLmr)
    if(is.null(errL)) {
      err <- paste0("ERR:02141:PCMBase:GaussianPCM.R:PCMLik:: There was a problem calculating the coefficients L,m,r. Error message from call to PCMLmr: ", PCMLmr, "; print(model):",
                    do.call(paste, c(as.list(capture.output(print(model))), list(sep="\n"))))
      errL <- PCMParseErrorMessage(err)
    } else {
      err <- PCMLmr
    }
    if(getOption("PCMBase.Errors.As.Warnings", TRUE)) {
      warning(err)
    } else {
      stop(err)
    }
    X0 <- model$X0
    attr(value.NA, 'X0') <- X0
    attr(value.NA, "error") <- errL

    return(value.NA)

  } else if(is.list(PCMLmr)) {

    L_root <- PCMLmr$L
    m_root <- PCMLmr$m
    r_root <- PCMLmr$r

    if(is.null(L_root) | is.null(m_root) | is.null(r_root)) {
      err <- paste0("ERR:02142:PCMBase:GaussianPCM.R:PCMLik:: The list returned by PCMLmr did not contain elements 'L', 'm' and 'r'.")
      if(getOption("PCMBase.Errors.As.Warnings", TRUE)) {
        warning(err)
      } else {
        stop(err)
      }

      errL <- PCMParseErrorMessage(err)

      X0 <- model$X0
      attr(value.NA, 'X0') <- X0
      attr(value.NA, "error") <- errL

      return(value.NA)

    }

    if(is.null(model$X0) || isTRUE(all(is.na(model$X0)))) {
      # set the root value to the one that maximizes the likelihood
      X0 <- try(solve(a=L_root + t(L_root), b = -m_root), silent = TRUE)
      if(class(X0) == "try-error") {
        err <- paste0(
          "ERR:02143:PCMBase:GaussianPCM.R:PCMLik:: There was a problem calculating X0 from the coefficients L,m,r. ", "L=", toString(L), "; m=", toString(m), "; r=", r,
          ". Error message from call to solve(a=L_root + t(L_root), b = -m_root):", X0, "; print(model):",
          do.call(paste, c(as.list(capture.output(print(model))), list(sep="\n"))))

        errL <- PCMParseErrorMessage(err)
        if(getOption("PCMBase.Errors.As.Warnings", TRUE)) {
          warning(err)
        } else {
          stop(err)
        }
        X0 <- NULL
        attr(value.NA, "X0") <- X0
        attr(value.NA, "error") <- errL

        return(value.NA)

      }

    } else {

      X0 <- model$X0

    }

    loglik <- try(X0 %*% L_root %*% X0 + m_root %*% X0 + r_root, silent = TRUE)
    if(class(loglik) == "try-error") {
      err <- paste0(
        "ERR:02144:PCMBase:GaussianPCM.R:PCMLik:: There was a problem calculating loglik from X0 and the coefficients L,m,r. ", "X0=", toString(X0), "L=", toString(L), "; m=", toString(m), "; r=", r,
        ". Error message from call to X0 %*% L_root %*% X0 + m_root %*% X0 + r_root:", loglik, "; print(model):",
        do.call(paste, c(as.list(capture.output(print(model))), list(sep="\n"))))

      errL <- PCMParseErrorMessage(err)
      if(getOption("PCMBase.Errors.As.Warnings", TRUE)) {
        warning(err)
      } else {
        stop(err)
      }

      attr(value.NA, "X0") <- X0
      attr(value.NA, "error") <- errL

      return(value.NA)

    }

    value <- try(as.vector(if(log) loglik else exp(loglik)), silent = TRUE)

    if(class(value) == "try-error") {
      err <- paste0(
        "ERR:02145:PCMBase:GaussianPCM.R:PCMLik:: There was a problem calculating value from loglik=", toString(loglik), ". Error message from call to as.vector(if(log) loglik else exp(loglik)):", value, "; print(model):",
        do.call(paste, c(as.list(capture.output(print(model))), list(sep="\n"))))

      errL <- PCMParseErrorMessage(err)
      if(getOption("PCMBase.Errors.As.Warnings", TRUE)) {
        warning(err)
      } else {
        stop(err)
      }

      attr(value.NA, "X0") <- X0
      attr(value.NA, "error") <- errL

      return(value.NA)

    } else if(is.na(value)) {

      err <- paste0(
        "ERR:02146:PCMBase:GaussianPCM.R:PCMLik:: There was a possible numerical problem, e.g. division of 0 by 0 when calculating the likelihood. value=", toString(value), "; calculated loglik=", toString(loglik), "; print(model):",
        do.call(paste, c(as.list(capture.output(print(model))), list(sep="\n"))), ". No error message was returned from the call to PCMLmr. Check for runtime warnings.")

      errL <- PCMParseErrorMessage(err)
      if(getOption("PCMBase.Errors.As.Warnings", TRUE)) {
        warning(err)
      } else {
        stop(err)
      }

      attr(value.NA, "X0") <- X0
      attr(value.NA, "error") <- errL

      return(value.NA)

    }

    # no errors detected, returning the calculated likelihood value:
    attr(value, "X0") <- X0
    return(value)
  }
}




#' Quadratic polynomial parameters A, b, C, d, E, f for each node
#' @description An S3 generic function that has to be implemented for every
#'  model class. This function is called by \code{\link{PCMLik}}.
#' @inheritParams PCMLik
#' @export
PCMAbCdEf <- function(tree, model, metaI=PCMInfo(NULL, tree, model, verbose), verbose = FALSE) {
  UseMethod("PCMAbCdEf", model)
}

#' @export
PCMAbCdEf.default <- function(tree, model, metaI=PCMInfo(NULL, tree, model, verbose), verbose = FALSE) {
  threshold_SV <- getOption("PCMBase.Threshold.SV", 1e-6)
  skip_singular <- getOption("PCMBase.Singular.Skip", TRUE)

  # number of regimes
  R <- metaI$RModel
  # number of tips
  N <- metaI$N
  # number of traits (variables)
  k <- metaI$k
  # number of nodes
  M <- metaI$M
  # present coordinates
  pc <- metaI$pc

  cond <- list()

  for(r in 1:R) {
    # create the conditional distribution for regime r
    cond[[r]] <- PCMCond(tree, model, r, metaI, verbose)
  }

  omega <- array(NA, dim=c(k, M))
  Phi <- array(NA, dim=c(k, k, M))
  V <- array(NA, dim=c(k, k, M))
  V_1 <- array(NA, dim=c(k, k, M))


  # returned general form parameters
  A <- array(NA, dim=c(k, k, M))
  b <- array(NA, dim=c(k, M))
  C <- array(NA, dim=c(k, k, M))
  d <- array(NA, dim=c(k, M))
  E <- array(NA, dim=c(k, k, M))
  f <- array(NA, dim=c(M))

  singular <- rep(FALSE, M)

  # vector of regime indices for each branch
  r <- metaI$r

  # identity k x k matrix
  I <- diag(k)

  # iterate over the edges
  for(edgeIndex in 1:(M-1)) {
    # parent node
    j <- tree$edge[edgeIndex, 1]
    # daughter node
    i <- tree$edge[edgeIndex, 2]

    # length of edge leading to i
    ti <- tree$edge.length[edgeIndex]

    # present coordinates in parent and daughte nodes
    kj <- pc[,j]
    ki <- pc[,i]

    omega[,i] <- cond[[r[edgeIndex]]]$omega(ti, edgeIndex, metaI)
    Phi[,,i] <- cond[[r[edgeIndex]]]$Phi(ti, edgeIndex, metaI)
    V[,,i] <- cond[[r[edgeIndex]]]$V(ti, edgeIndex, metaI)


    # check that V[ki,ki,] is non-singular
    svdV = svd(matrix(V[ki,ki,i], sum(ki)), 0, 0)$d

    if(is.na(min(svdV)/max(svdV)) || min(svdV)/max(svdV) < threshold_SV) {
      singular[edgeIndex] <- TRUE
      if(!skip_singular) {
        err <- paste0(
          "ERR:02131:PCMBase:GaussianPCM.R:PCMAbCdEf.default:",i,":",
          " The matrix V for node ", i,
          " is nearly singular: min(svdV)/max(svdV)=", min(svdV)/max(svdV),
          ". Check the model parameters and the length of the branch",
          " leading to the node. For details on this error, read the User Guide.")
        stop(err)
      }
    }

    if(!singular[edgeIndex]) {
      V_1[ki,ki,i] <- solve(V[ki,ki,i])

      `%op%` <- if(sum(ki) > 1) `%*%` else `*`

      A[ki,ki,i] <- (-0.5*V_1[ki,ki,i])
      E[kj,ki,i] <- t(Phi[ki,kj,i]) %op% V_1[ki,ki,i]
      b[ki,i] <- V_1[ki,ki,i] %*% omega[ki,i]
      C[kj,kj,i] <- -0.5 * matrix(E[kj,ki,i], sum(kj), sum(ki)) %*% matrix(Phi[ki,kj,i], sum(ki), sum(kj))
      d[kj,i] <- -E[kj,ki,i] %op% omega[ki,i]
      f[i] <- -0.5 * (t(omega[ki,i]) %*% V_1[ki,ki,i] %*% omega[ki,i] + sum(ki)*log(2*pi) + log(det(as.matrix(V[ki,ki,i]))))
    }
  }

  list(A=A, b=b, C=C, d=d, E=E, f=f, omega=omega, Phi=Phi, V=V, V_1=V_1, singular = singular)
}

#' Quadratic polynomial parameters L, m, r
#'
#' @description
#'
#' @inheritParams PCMLik
#' @param root.only logical indicatin whether to return the calculated values of L,m,r
#'  only for the root or for all nodes in the tree.
#' @return A list with the members A,b,C,d,E,f,L,m,r for all nodes in the tree or
#'   only for the root if root.only=TRUE.
#' @export
PCMLmr <- function(X, tree, model,
                metaI = PCMInfo(X, tree, model, verbose = verbose),
                root.only = TRUE, verbose = FALSE) {
  UseMethod("PCMLmr", metaI)
}


#' @export
PCMLmr.default <- function(
  X, tree, model,
  metaI = PCMInfo(X, tree, model, verbose = verbose),
  root.only = FALSE,
  verbose = FALSE) {

  unJ <- 1

  N <- metaI$N; M <- metaI$M; k <- metaI$k;

  edge <- tree$edge
  endingAt <- metaI$endingAt
  nodesVector <- metaI$nodesVector
  nodesIndex <- metaI$nodesIndex
  nLevels <- metaI$nLevels
  unVector <- metaI$unVector
  unIndex <- metaI$unIndex
  pc <- metaI$pc

  L <- array(0, dim=c(k, k, M))
  m <- array(0, dim=c(k, M))
  r <- array(0, dim=c(M))


  AbCdEf <- PCMAbCdEf(tree = tree, model = model, metaI = metaI, verbose = verbose)


  # avoid redundant calculation
  log2pi <- log(2*pi)

  for(level in 1:nLevels) {
    nodes <- nodesVector[(nodesIndex[level]+1):nodesIndex[level+1]]
    es <- endingAt[nodes]

    if(nodes[1] <= N) {
      # all es pointing to tips

      for(edgeIndex in es) {
        if(!AbCdEf$singular[edgeIndex]) {
          # parent and daughter nodes
          j <- edge[edgeIndex, 1]; i <- edge[edgeIndex, 2];
          # present coordinates
          kj <- pc[, j]; ki <- pc[, i];

          # ensure symmetry of L[,,i]
          L[,,i] <- 0.5 * (AbCdEf$C[,,i] + t(AbCdEf$C[,,i]))

          r[i] <- with(AbCdEf, t(X[ki,i]) %*% A[ki,ki,i] %*% X[ki,i] +
                         t(X[ki,i]) %*% b[ki,i] + f[i])

          m[kj,i] <- with(AbCdEf, d[kj,i] + matrix(E[kj,ki,i], sum(kj), sum(ki)) %*% X[ki,i])
        }
      }
    } else {
      # edges pointing to internal nodes, for which all children
      # nodes have been visited
      for(edgeIndex in es) {
        if(!AbCdEf$singular[edgeIndex]) {
          # parent and daughter nodes
          j <- edge[edgeIndex, 1]; i <- edge[edgeIndex, 2];
          # present coordinates
          kj <- pc[, j]; ki <- pc[, i];

          # auxilary variables to avoid redundant evaluation
          AplusL <- as.matrix(AbCdEf$A[ki,ki,i] + L[ki,ki,i])
          # ensure symmetry of AplusL, this should guarantee that AplusL_1 is symmetric
          # as well (unless solve-implementation is buggy.)
          AplusL <- 0.5 * (AplusL + t(AplusL))

          AplusL_1 <- solve(AplusL)

          EAplusL_1 <- matrix(AbCdEf$E[kj,ki,i], sum(kj), sum(ki)) %*% AplusL_1
          logDetVNode <- log(det(-2*AplusL))

          # here it is important that we first evaluate r[i] and then m[i,kj]
          # since the expression for r[i] refers to to the value of m[i,ki]
          # before updating it.
          r[i] <- with(AbCdEf, f[i]+r[i]+(sum(ki)/2)*log2pi-.5*logDetVNode -
                         .25*t(b[ki,i]+m[ki,i]) %*% AplusL_1 %*% (b[ki,i]+m[ki,i]))

          m[kj,i] <- with(AbCdEf, d[kj,i] - .5*EAplusL_1 %*% (b[ki,i]+m[ki,i]))

          L[kj,kj,i] <- with(
            AbCdEf,
            C[kj,kj,i] -.25*EAplusL_1 %*% t(matrix(E[kj,ki,i], sum(kj), sum(ki))))

          # ensure symmetry of L:
          L[kj,kj,i] <- 0.5 * (L[kj,kj,i] + t(L[kj,kj,i]))
        }
      }
    }

    # add up to parents
    while(length(es)>0) {
      un <- unVector[(unIndex[unJ]+1):unIndex[unJ+1]]
      unJ <- unJ+1
      L[,,edge[es[un], 1]] <- L[,,edge[es[un], 1]] + L[,,edge[es[un], 2]]
      m[,edge[es[un], 1]] <- m[,edge[es[un], 1]] + m[,edge[es[un], 2]]
      r[edge[es[un], 1]] <- r[edge[es[un], 1]] + r[edge[es[un], 2]]
      es <- es[-un]
    }
  }

  if(root.only) {
    list(L = L[,,N+1],
         m = m[,N+1],
         r = r[N+1])
  } else {
    c(AbCdEf[c("A", "b", "C", "d", "E", "f", "V", "V_1")],
      list(L = L, m = m, r = r))
  }
}

#' Sums of pairs of elements in a vector
#' @param lambda a numeric vector
#' @return a squared symmetric matrix with elem_ij=lambda_i+lambda_j.
#'
PCMPairSums <- function(lambda) {
  sapply(lambda, function(li) sapply(lambda, function(lj) li+lj))
}

#' Eigen-decomposition of a matrix H
#' @param H a numeric matrix
#' @return a list with elements as follows:
#' \item{lambda}{a vector of the eigenvalues of H}
#' \item{P}{a squared matrix with column vectors, the eigenvectors of H corresponding to the
#' eigenvalues in \code{lambda}}
#' \item{P_1}{the inverse matrix of \code{P}}.
#' @details The function fails with an error message if H is defective, that is, if its matrix of
#' eigenvectors is computationally singular. The test for singularity is based on the \code{\link{rcond}} function.
#'
PCMPLambdaP_1 <- function(H) {
  # here the argument H is a matrix specifying the alphas in a OU process
  r <- eigen(H)
  if(isTRUE(all.equal(rcond(r$vectors),0))) {
    stop(paste0("ERR:02141:PCMBase:GaussianPCM.R:PCMPLambdaP_1:: The provided H matrix is defective. Its matrix of eigenvectors is computationally singular. H=", toString(H)))
  }
  list(lambda=r$values, P=r$vectors, P_1=solve(r$vectors))
}


#' Create a function of time that calculates (1-exp(-lambda_{ij}*time))/lambda_{ij}
#' for every element lambda_{ij} of the input matrix Lambda_ij.
#'
#' @param Lambda_ij a squared numerical matrix of dimension k x k
#' @param threshold.Lambda_ij a 0-threshold for abs(Lambda_i + Lambda_j), where Lambda_i
#'   and Lambda_j are eigenvalues of the parameter matrix H. This
#'   threshold-value is used as a condition to
#'   take the limit time of the expression `(1-exp(-Lambda_ij*time))/Lambda_ij` as
#'   `(Lambda_i+Lambda_j) --> 0`. You can control this value by the global option
#'   "PCMBase.Threshold.Lambda_ij". The default value (1e-8) is suitable for branch
#'    lengths bigger than 1e-6. For smaller branch lengths, you may want to
#'    increase the threshold value using, e.g.
#'   `options(PCMBase.Threshold.Lambda_ij=1e-6)`.
#' @details the function (1-exp(-lambda_{ij}*time))/lambda_{ij} corresponds to the product
#' of the CDF of an exponential distribution with rate Lambda_{ij} multiplied by its mean value (mean waiting time).
#' @return a function of time returning a matrix with entries formed from the
#'  above function or the limit, time, if |Lambda_{ij}|<=trehshold0.
PCMPExpxMeanExp <- function(
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


#' Variance-covariance matrix of an OU process with optional measurement error and jump at the start
#' @param H a numerical k x k matrix - selection strength parameter.
#' @param Sigma a numerical k x k matrix - neutral dift unit-time variance-covariance matrix.
#' @param Sigmaj is the variance matrix of the normal jump distribution (default is NULL).
#' @param xi a vector of 0's and 1's corresponding to each branch in the tree. A value of 1
#' indicates that a jump takes place at the beginning of the branch. This arugment is only
#' used if Sigmaj is not NULL. Default is NULL.
#'
#' @param threshold.Lambda_ij a 0-threshold for abs(Lambda_i + Lambda_j), where Lambda_i
#' and Lambda_j are eigenvalues of the parameter matrix H. This threshold-values is used as
#' a condition to take the limit time of the expression `(1-exp(-Lambda_ij*time))/Lambda_ij`
#' as `(Lambda_i+Lambda_j) --> 0`. You can control this value by the global option
#' "PCMBase.Threshold.Lambda_ij". The default value (1e-8) is suitable for branch lengths
#' bigger than 1e-6. For smaller branch lengths, you may want to increase the threshold
#' value using, e.g.  `options(PCMBase.Threshold.Lambda_ij=1e-6)`.
#' @return a function of one numerical argument (time) and an integer indicating the branch-index
#' that is used to check the corresponding element in xi.
#' @export
PCMCondVOU <- function(
  H, Sigma , Sigmae = NULL, Sigmaj = NULL, xi = NULL,
  e_Ht = NULL,
  threshold.Lambda_ij = getOption("PCMBase.Threshold.Lambda_ij", 1e-8)) {

  force(H)
  force(Sigma)
  force(Sigmae)
  force(Sigmaj)
  force(xi)
  force(e_Ht)

  if(is.null(dim(H)) | is.null(dim(Sigma))) {
    stop('ERR:02151:PCMBase:GaussianPCM.R:PCMCondVOU:: H and Sigma must be k x k matrices.')
  }

  k <- dim(Sigma)[1]

  if(!is.matrix(Sigma) | !is.matrix(H) | !isTRUE(all.equal(c(dim(H), dim(Sigma)), rep(k, 4)))) {
    stop(paste0("ERR:02152:PCMBase:GaussianPCM.R:PCMCondVOU:: H and Sigma must be  ", k, " x ", k, " matrices."))
  }

  if(!is.null(Sigmaj) & !is.matrix(Sigmaj) & !isTRUE(all.equal(dim(Sigmaj), rep(k, 2)))) {
    stop(paste0("ERR:02153:PCMBase:GaussianPCM.R:PCMCondVOU:: Sigmaj must be NULL or a ", k, " x ", k, " matrix."))
  }

  if(!is.null(e_Ht) & !is.matrix(e_Ht) & !isTRUE(all.equal(dim(e_Ht), rep(k, 2)))) {
    stop(paste0("ERR:02154:PCMBase:GaussianPCM.R:PCMCondVOU:: e_Ht must be NULL or a ", k, " x ", k, " matrix."))
  }

  PLP_1 <- PCMPLambdaP_1(H)

  Lambda_ij <- PCMPairSums(PLP_1$lambda)
  fLambda_ij <- PCMPExpxMeanExp(Lambda_ij, threshold.Lambda_ij)
  P_1SigmaP_t <- PLP_1$P_1 %*% Sigma %*% t(PLP_1$P_1)


  function(t, edgeIndex, metaI, e_Ht = NULL) {
    res <- PLP_1$P %*% (fLambda_ij(t) * P_1SigmaP_t) %*% t(PLP_1$P)
    if(!is.null(Sigmaj)) {
      if(is.null(e_Ht)) {
        e_Ht <- expm(-t*H)
      }
      res <- res + xi[edgeIndex]*(e_Ht %*% Sigmaj %*% t(e_Ht))
    }
    if(!is.null(Sigmae) & metaI$edge[edgeIndex,2] <= metaI$N) {
      res <- res + Sigmae
    }
    Re(res)
  }
}

PCMCondRandom <- function(PCMCondObject, n=1, x0, t, edgeIndex, metaI) {
  with(PCMCondObject, {
    rmvnorm(n=n, mean = omega(t, edgeIndex, metaI) + Phi(t, edgeIndex, metaI)%*%x0, sigma=V(t, edgeIndex, metaI))
  })
}

PCMCondDensity <- function(PCMCondObject, x, x0, t, edgeIndex, metaI, log=FALSE) {
  with(PCMCondObject, {
    dmvnorm(x, mean=omega(t, edgeIndex, metaI) + Phi(t, edgeIndex, metaI)%*%x0, sigma=V(t, edgeIndex, metaI), log=log)
  })
}
