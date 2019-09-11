# Copyright 2016-2019 Venelin Mitov
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
PCMCond.GaussianPCM <- function(
  tree, model, r=1,
  metaI = PCMInfo(NULL, tree, model, verbose = verbose), verbose = FALSE) {
  stop(paste(
    'ERR:02111:PCMBase:GaussianPCM.R:PCMCond.GaussianPCM:: This is an abstract',      'function which only defines an interface but should not be called ',
    'explicitly. Possibly you forgot implementing PCMCond for a daughter class.'
    ))
}

#' @export
PCMMean.GaussianPCM <- function(
  tree, model, X0 = model$X0,
  metaI = PCMInfo(
    X = NULL, tree = tree, model = model, verbose = verbose),
  internal = FALSE, verbose = FALSE)  {
  if(is.Transformable(model)) {
    model <- PCMApplyTransformation(model)
  }
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

  # preorder order of the edges in the tree
  preord <- metaI$preorder

  # create a list of random generator functions for each regime
  PCMCondObjects <- lapply(seq_len(metaI$RModel), function(r) {

      PCMCond(tree, model = model, r = r, metaI = metaI, verbose = verbose)
  })

  Mu <- matrix(as.double(NA), k, M)
  Mu[, N+1] <- X0

  for(edgeIndex in preord) {
    cond <- PCMCondObjects[[metaI$r[edgeIndex]]]
    t <- tree$edge.length[edgeIndex]
    # parent node
    j <- tree$edge[edgeIndex, 1L]
    # daughter node
    i <- tree$edge[edgeIndex, 2L]
    Mu[, i] <- cond$omega(t, edgeIndex, metaI) + cond$Phi(t, edgeIndex, metaI) %*% Mu[, j]
  }

  if(internal) {
    Mu
  } else {
    Mu[, 1:N]
  }
}

#' @importFrom ape mrca
#' @export
PCMVar.GaussianPCM <- function(
  tree, model,
  W0 = matrix(0.0, PCMNumTraits(model), PCMNumTraits(model)),
  SE = matrix(0.0, PCMNumTraits(model), PCMTreeNumTips(tree)),
  metaI = PCMInfo(
    X = NULL, tree = tree, model = model, SE = SE, verbose = verbose),
  internal = FALSE, verbose = FALSE)  {

  threshold_SV <- metaI$PCMBase.Threshold.SV
  skip_singular <- metaI$PCMBase.Skip.Singular
  threshold_skip_singular <- metaI$PCMBase.Threshold.Skip.Singular
  threshold_eigval <- metaI$PCMBase.Threshold.EV

  if(is.Transformable(model)) {
    model <- PCMApplyTransformation(model)
  }
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

  # preorder order of the edges in the tree
  preord <- metaI$preorder

  # create a list of random generator functions for each regime
  PCMCondObjects <- lapply(1:metaI$RModel, function(r) {
    PCMCond(tree, model = model, r = r, metaI = metaI, verbose = verbose)
  })

  BlockI <- function(i) {
    (i-1)*k + seq_len(k)
  }

  # Variance-covariance matrix at the tips
  W <- matrix(as.double(NA), k*N, k*N)

  # Variance-blocks matrix for all nodes
  Wii <- matrix(as.double(NA), k, k*M)
  Wii[, BlockI(N+1)] <- W0

  # need to set the names of the nodes and tips to their integer indices
  # converted to character strings this affects the local tree
  PCMTreeSetLabels(tree)
  listPathsToRoot <- PCMTreeListRootPaths(tree)

  MRCA <- function(i, j) {
    if(i == j) {
      i
    } else {
      path_i <- listPathsToRoot[[i]]
      path_j <- listPathsToRoot[[j]]

      ii <- length(path_i)
      jj <- length(path_j)

      # the root is certainly a common ancestor
      a <- as.integer(NA)
      while(ii > 0 && jj > 0) {
        if(path_i[[ii]] != path_j[[jj]]) {
          break
        } else {
          a <- path_i[[ii]]
          ii <- ii - 1
          jj <- jj - 1
        }
      }
      a
    }
  }

  # Products of the Phi matrix in the direction from the nodes to the root
  # listProdPhi[[i]] = Phi_i %*% Phi_{parent(i)} %*% ... %*% I(root)
  listProdPhi <- list()

  # first we calculate the variance (diagonal) blocks (i,i) following preorder traversal
  for(edgeIndex in preord) {
    cond <- PCMCondObjects[[metaI$r[edgeIndex]]]
    t <- tree$edge.length[edgeIndex]
    # parent node
    j <- tree$edge[edgeIndex, 1L]
    # daughter node
    i <- tree$edge[edgeIndex, 2L]

    Phi_i <- cond$Phi(t, edgeIndex, metaI)

    # handle measurement error
    if(i <= metaI$N) {
      # tip node
      VE <- metaI$VE[, , i]
    } else {
      # internal node
      VE <- matrix(0.0, metaI$k, metaI$k)
    }

    V <- cond$V(t, edgeIndex, metaI) + VE

    # force symmetric V
    V <- 0.5 * (V + t(V))

    # check that V is non-singular
    svdV = svd(V, 0, 0)$d

    eigval <- eigen(V, symmetric = TRUE, only.values = TRUE)$values

    if(is.na(min(svdV)/max(svdV)) ||
       min(svdV)/max(svdV) < threshold_SV ||
       eigval[k] < threshold_eigval) {
      if( !skip_singular ||
          # We try to support singular branches leading to tips, but we don't
          # support such during a likelihood calculation (Function PCMAbCdEf).
          # i <= metaI$N ||
          t > threshold_skip_singular ) {
        err <- paste0(
          "PCMVar.GaussianPCM:",i,":",
          " The matrix V for node ", i, " (branch length=", t, ")",
          " is nearly singular: min(svdV)/max(svdV)=", min(svdV)/max(svdV),
          ". smallest eigenvalue: ", eigval[k],
          "; Check PCMOptions()$PCMBase.Threshold.EV, ", "
          PCMOptions()$PCMBase.Threshold.SV, ",
          "PCMOptions()$PCMBase.Threshold.Skip.Singular ",
          "and the model parameters and the length of the branch ",
          "leading to the node.")
        stop(err)
      } else {
        # skip the singular branch by taking the values from j unchanged
        Wii[, BlockI(i)] <- Wii[, BlockI(j)]

        listProdPhi[[i]] <- lapply(listPathsToRoot[[i]], function(jAnc) {
          if(jAnc == j) {
            diag(1, k, k)
          } else {
            listProdPhi[[j]][[as.character(jAnc)]]
          }
        })
      }
    } else {
      Wii[, BlockI(i)] <- Phi_i %*% Wii[, BlockI(j)] %*% t(Phi_i) + V

      # force symmetry for Wii[, BlockI(i)]
      Wii[, BlockI(i)] <- 0.5*(Wii[, BlockI(i)] + t(Wii[, BlockI(i)]))

      listProdPhi[[i]] <- lapply(listPathsToRoot[[i]], function(jAnc) {
        if(jAnc == j) {
          Phi_i
        } else {
          Phi_i %*% listProdPhi[[j]][[as.character(jAnc)]]
        }
      })
    }
  }

  if(N > 1) {
    for(i in seq_len(N)) {
      for(j in seq_len(N)) {
        mrca_ij <- MRCA(i, j)

        if(mrca_ij == i) {
          ProdPhi_i <- diag(k)
        } else {
          ProdPhi_i <- listProdPhi[[i]][[as.character(mrca_ij)]]
        }

        if(mrca_ij == j) {
          ProdPhi_j <- diag(k)
        } else {
          ProdPhi_j <- listProdPhi[[j]][[as.character(mrca_ij)]]
        }

        W[BlockI(i), BlockI(j)] <- # W[BlockI(j), BlockI(i)] <-
           ProdPhi_i %*% Wii[, BlockI(mrca_ij)] %*% t(ProdPhi_j)

        # assign mirror block and
        # force symmetry for W[BlockI(i), BlockI(j)] and W[BlockI(j), BlockI(i)]
        W[BlockI(i), BlockI(j)] <- W[BlockI(j), BlockI(i)] <-
          0.5 * (W[BlockI(i), BlockI(j)] + t(W[BlockI(i), BlockI(j)]))
      }
    }
  } else if(N == 1) {
    W[BlockI(1), BlockI(1)] <- Wii[, BlockI(1)]
  }

  if(internal) {
    #list(W = 0.5*(W + t(W)), Wii = Wii)
    list(W = W, Wii = Wii)
  } else {
    #0.5*(W + t(W))
    W
  }
}


#' @export
PCMSim.GaussianPCM <- function(
  tree, model, X0,
  SE = matrix(0.0, PCMNumTraits(model), PCMTreeNumTips(tree)),
  metaI = PCMInfo(
    X = NULL, tree = tree, model = model, SE = SE, verbose = verbose),
  verbose = FALSE) {

  if(is.Transformable(model)) {
    model <- PCMApplyTransformation(model)
  }

  if(length(X0)!=metaI$k) {
    stop(paste('GaussianPCM.R:PCMSim:: X0 must be of length', metaI$k, '.'))
  }

  values <- matrix(0, nrow=metaI$k, ncol=dim(tree$edge)[1]+1)
  values[, metaI$N + 1] <- X0

  preord <- metaI$preorder

  # create a list of random generator functions for each regime
  PCMCondObjects <- lapply(1:metaI$RModel, function(r) {
    PCMCond(tree, model = model, r = r, metaI = metaI, verbose = verbose)
  })

  for(edgeIndex in preord) {
    obj <- PCMCondObjects[[metaI$r[edgeIndex]]]
    # parent node
    j <- tree$edge[edgeIndex, 1]
    # daughter node
    i <- tree$edge[edgeIndex, 2]

    if(i <= metaI$N) {
      # daughter is a tip
      VE <- metaI$VE[, , i]
    } else {
      # daughter is an internal node
      VE <- matrix(0.0, metaI$k, metaI$k)
    }

    if(!is.null(obj$random)) {
      values[, i] <-
        obj$random(
          n=1,
          x0 = values[, j],
          t = tree$edge.length[edgeIndex],
          edgeIndex = edgeIndex,
          VE = VE)
    } else {
      values[, i] <-
        PCMCondRandom(
          obj,
          n=1,
          x0 = values[, j],
          t = tree$edge.length[edgeIndex],
          edgeIndex = edgeIndex,
          metaI = metaI,
          VE = VE)
    }
  }

  values
}

#' @importFrom utils capture.output
#'
#' @export
PCMLik.GaussianPCM <- function(
  X, tree, model,
  SE = matrix(0.0, PCMNumTraits(model), PCMTreeNumTips(tree)),
  metaI = PCMInfo(
    X = X, tree = tree, model = model, SE = SE, verbose = verbose),
  log = TRUE,
  verbose = FALSE) {

  if(is.Transformable(model)) {
    model <- PCMApplyTransformation(model)
  }

  if(is.function(metaI)) {
    metaI <- metaI(
      X = X, tree = tree, model = model, SE = SE, verbose = verbose)
  }

  # will change this value if there is no error
  value.NA <- getOption("PCMBase.Value.NA", as.double(NA))

  PCMLmr <- try(PCMLmr(
    X = X, tree = tree, model = model, SE = SE, metaI = metaI,
    verbose = verbose, root.only = TRUE),
    silent = TRUE)

  if(class(PCMLmr) == "try-error") {
    errL <- PCMParseErrorMessage(PCMLmr)
    if(is.null(errL)) {
      err <- paste0("GaussianPCM.R:PCMLik:: There was a problem calculating the coefficients L,m,r. Error message from call to PCMLmr: ", PCMLmr, "; print(model):",
                    do.call(paste, c(as.list(capture.output(print(model))), list(sep="\n"))))
      errL <- PCMParseErrorMessage(err)
    } else {
      err <- PCMLmr
    }

    if(getOption("PCMBase.Raise.Lik.Errors", TRUE)) {
      if(getOption("PCMBase.Errors.As.Warnings", TRUE)) {
        warning(err)
      } else {
        stop(err)
      }
    }


    k0 <- metaI$pc[, metaI$N + 1L]
    X0 <- rep(NaN, metaI$k)
    X0[k0] <- model$X0[k0]

    attr(value.NA, 'X0') <- X0
    attr(value.NA, "error") <- errL

    return(value.NA)

  } else if(is.list(PCMLmr)) {
    L0 <- as.matrix(PCMLmr$L)
    m0 <- PCMLmr$m
    r0 <- PCMLmr$r
    k0 <- metaI$pc[, metaI$N + 1L]

    if(is.null(L0) | is.null(m0) | is.null(r0)) {
      err <- paste0("GaussianPCM.R:PCMLik:: The list returned by PCMLmr did not contain elements 'L', 'm' and 'r'.")

      if(!getOption("PCMBase.Ignore.Lmr.Errors", FALSE)) {
        if(getOption("PCMBase.Errors.As.Warnings", TRUE)) {
          warning(err)
        } else {
          stop(err)
        }
      }

      errL <- PCMParseErrorMessage(err)

      k0 <- metaI$pc[, metaI$N + 1L]
      X0 <- rep(NaN, metaI$k)
      X0[k0] <- model$X0[k0]

      attr(value.NA, 'X0') <- X0
      attr(value.NA, "error") <- errL

      return(value.NA)
    }

    if(is.null(model$X0) || isTRUE(all(is.na(model$X0)))) {
      # set the root value to the one that maximizes the likelihood
      X0 <- rep(NaN, metaI$k)
      X0[k0] <- rep(NA_real_, sum(k0))

      X0[k0] <- try(solve(
        a=L0[k0,k0,drop=FALSE] + t(L0[k0,k0,drop=FALSE]),
        b = -m0[k0]), silent = TRUE)
      if(class(X0[[1]]) == "try-error") {
        err <- paste0(
          "GaussianPCM.R:PCMLik:: There was a problem calculating X0 from L0,m0,r0. ",
          "L0=", toString(L0), "; m0=", toString(m0),
          "; r0=", r0, ". Error message:", X0[[1]], "\n")

        errL <- PCMParseErrorMessage(err)
        if(getOption("PCMBase.Raise.Lik.Errors", TRUE)) {
          if(getOption("PCMBase.Errors.As.Warnings", TRUE)) {
            warning(err)
          } else {
            stop(err)
          }
        }
        X0 <- NULL
        attr(value.NA, "X0") <- X0
        attr(value.NA, "error") <- errL

        return(value.NA)

      }

    } else {
      k0 <- metaI$pc[, metaI$N + 1L]
      X0 <- rep(NaN, metaI$k)
      X0[k0] <- model$X0[k0]
    }

    loglik <- try(X0[k0] %*% L0[k0,k0,drop=FALSE] %*% X0[k0] + m0[k0] %*% X0[k0] + r0, silent = TRUE)
    if(class(loglik) == "try-error") {
      err <- paste0(
        "GaussianPCM.R:PCMLik:: There was a problem calculating loglik from X0 and the coefficients L,m,r. ", "X0=", toString(X0), "L0=", toString(L0), "; m0=", toString(m0), "; r0=", r0,
        ". Error message from call to X0 %*% L0 %*% X0 + m0 %*% X0 + r0:", loglik, "\n")

      errL <- PCMParseErrorMessage(err)
      if(getOption("PCMBase.Raise.Lik.Errors", TRUE)) {
        if(getOption("PCMBase.Errors.As.Warnings", TRUE)) {
          warning(err)
        } else {
          stop(err)
        }
      }

      attr(value.NA, "X0") <- X0
      attr(value.NA, "error") <- errL

      return(value.NA)

    }

    value <- try(as.vector(if(log) loglik else exp(loglik)), silent = TRUE)

    if(class(value) == "try-error") {
      err <- paste0(
        "GaussianPCM.R:PCMLik:: There was a problem calculating value from loglik=", toString(loglik), ". Error message from call to as.vector(if(log) loglik else exp(loglik)):", value, "; print(model):",
        do.call(paste, c(as.list(capture.output(print(model))), list(sep="\n"))))

      errL <- PCMParseErrorMessage(err)

      if(getOption("PCMBase.Raise.Lik.Errors", TRUE)) {
        if(getOption("PCMBase.Errors.As.Warnings", TRUE)) {
          warning(err)
        } else {
          stop(err)
        }
      }

      attr(value.NA, "X0") <- X0
      attr(value.NA, "error") <- errL

      return(value.NA)

    } else if(is.na(value)) {
      err <- paste0(
        "GaussianPCM.R:PCMLik:: There was a possible numerical problem, e.g. division of 0 by 0 when calculating the likelihood. value=", toString(value), "; calculated loglik=", toString(loglik), "; print(model):",
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

#' @importFrom data.table setnames
#' @importFrom utils tail
#' @export
PCMLikTrace.GaussianPCM <- function(
  X, tree, model,
  SE = matrix(0.0, PCMNumTraits(model), PCMTreeNumTips(tree)),
  metaI = PCMInfo(
    X = X, tree = tree, model = model, SE = SE, verbose = verbose),
  log = TRUE,
  verbose = FALSE
) {
  if(is.Transformable(model)) {
    model <- PCMApplyTransformation(model)
  }

  if(is.function(metaI)) {
    metaI <- metaI(
      X = X, tree = tree, model = model, SE = SE, verbose = verbose)
  }

  trace <- PCMLmr(X = X, tree, model, metaI = metaI)
  names(trace) <- sapply(names(trace), function(name) {
    if(name %in% c("L", "m", "r")) {
      paste0(name, "_{ji}")
    } else {
      paste0(name, "_i")
    }
  })

  traceList <- lapply(
    names(trace), function(name) {
      res <- if(is.list(trace[[name]])) {
        trace[[name]]
      } else if(is.array(trace[[name]]) && length(dim(trace[[name]])) == 3) {
        lapply(seq_len(dim(trace[[name]])[3]), function(i) {
          trace[[name]][,,i]
        })
      } else if(is.array(trace[[name]]) && length(dim(trace[[name]])) == 2) {
        lapply(seq_len(dim(trace[[name]])[2]), function(i) {
          trace[[name]][,i]
        })
      } else if(
        is.array(trace[[name]]) && length(dim(trace[[name]])) == 1 ||
        is.vector(trace[[name]])) {

        trace[[name]]
      }
    })
  names(traceList) <- names(trace)

  traceTable <- do.call(data.table, traceList)

  # prevent warning during R CMD CHECK
  j <- .I <- i <- t_i <- k_i <- X_i <- .N <- VE_i <-
    omega_i <- Phi_i <- V_i <- L_i <- m_i <-
    r_i <- nodeId <- `L_{ji}` <- `m_{ji}` <- `r_{ji}` <- ff_i <- f_i <-
    `rr_{ji}` <- L_i <- m_i <- r_i <- `\\hat{X}_i` <- `\\ell\\ell_i` <- NULL

  traceTable[, j:=sapply(.I, function(nodeId) {
    if(nodeId != metaI$N+1) {
      PCMTreeGetLabels(tree)[PCMTreeGetParent(tree, nodeId)]
    } else {
      "ND"
    }
  })]
  traceTable[, i:=PCMTreeGetLabels(tree)[.I]]
  traceTable[, t_i:=lapply(.I, function(.) tree$edge.length[match(., tree$edge[, 2])])]
  traceTable[metaI$N+1]$t_i[[1]] <- list("ND")

  traceTable[, k_i:=lapply(
    seq_len(ncol(metaI$pc)), function(i) which(metaI$pc[,i]))]

  traceTable[, X_i:=lapply(seq_len(.N), function(nodeId) {
    nodeLabel <- i[nodeId]
    x <- if(nodeId <= metaI$N) {
      X[, nodeId]
    } else if(nodeId == metaI$N+1) {
      x <- rep(NaN, PCMNumTraits(model))
      if(any(is.finite(k_i[[nodeId]]))) {
        x[k_i[[nodeId]]] <- model$X0[k_i[[nodeId]]]
      }
      x
    } else {
      x <- rep(NaN, PCMNumTraits(model))
      if(any(is.finite(k_i[[nodeId]]))) {
        x[k_i[[nodeId]]] <- NA_real_
      }
      x
    }
    x <- as.character(x)
    mat <- xtable(
      as.matrix(x), align=rep("", ncol(as.matrix(x)) + 1),
      digits = getOption("digits", default = 2))
    latex <- capture.output(
      print(mat, floating=FALSE, tabular.environment="bmatrix",
            hline.after=NULL, include.rownames=FALSE,
            include.colnames=FALSE, comment = FALSE,
            NA.string = "NA"))
    paste0(" $", do.call(paste0, as.list(latex)), "$ ")
  })]

  traceTable[, VE_i:=lapply(seq_len(.N), function(nodeId) {
    nodeLabel <- i[nodeId]
    ve <- matrix(NA_real_, metaI$k, metaI$k)
    ve[k_i[[nodeId]], k_i[[nodeId]]] <- if(nodeId <= metaI$N) {
      metaI$VE[k_i[[nodeId]], k_i[[nodeId]], nodeId]
    } else {
      matrix(0, length(k_i[[nodeId]]), length(k_i[[nodeId]]))
    }
    ve
  })]

  setcolorder(traceTable, neworder = tail(seq_len(ncol(traceTable)), n = 5))

  traceTable[, omega_i:=lapply(.I, function(nodeId) {
    if(nodeId == metaI$N+1) {
      "ND"
    } else {
      vec <- rep(NA_real_, metaI$k)
      vec[k_i[[nodeId]]] <- omega_i[[nodeId]][k_i[[nodeId]]]
      vec
    }
  })]
  traceTable[, Phi_i:=lapply(.I, function(nodeId) {
    if(nodeId == metaI$N+1) {
      "ND"
    } else {
      k_j <- k_i[[PCMTreeGetParent(tree, nodeId)]]
      mat <- matrix(NA_real_, metaI$k, metaI$k)
      mat[k_i[[nodeId]], k_j] <- Phi_i[[nodeId]][k_i[[nodeId]], k_j]
      mat
    }
  })]
  traceTable[, V_i:=lapply(.I, function(nodeId) {
    if(nodeId == metaI$N+1) {
      "ND"
    } else {
      mat <- matrix(NA_real_, metaI$k, metaI$k)
      mat[k_i[[nodeId]], k_i[[nodeId]]] <-
        V_i[[nodeId]][k_i[[nodeId]], k_i[[nodeId]]]
      mat
    }
  })]
  traceTable[, L_i:=lapply(.I, function(nodeId) {
    if(nodeId <= metaI$N) {
      "ND"
    } else {
      mat <- matrix(NA_real_, metaI$k, metaI$k)
      mat[k_i[[nodeId]], k_i[[nodeId]]] <- 0.0
      mat
    }
  })]
  traceTable[, m_i:=lapply(.I, function(nodeId) {
    if(nodeId <= metaI$N) {
      "ND"
    } else {
      vec <- rep(NA_real_, metaI$k)
      vec[k_i[[nodeId]]] <- 0.0
      vec
    }
  })]
  traceTable[, r_i:=lapply(.I, function(nodeId) {
    if(nodeId <= metaI$N) {
      "ND"
    } else {
      0.0
    }
  })]

  # nodeId column needed because setting a key will sort the table on the
  # key-column.
  traceTable[, nodeId:=.I]
  setkey(traceTable, i)
  for(rowId in seq_len(nrow(traceTable))) {
    # daughter node label
    ii <- traceTable[rowId, i]
    # parent node label
    jj <- traceTable[rowId, j]

    if(traceTable[rowId, nodeId] != metaI$N + 1) {
      # not at the root node, so a parent node certainly exists
      traceTable[
        list(jj),
        L_i:=list(list( as.matrix(L_i[[1]] + traceTable[list(ii), `L_{ji}`[[1]]]) ))]
      traceTable[
        list(jj),
        m_i:=list(list(m_i[[1]] + traceTable[list(ii), `m_{ji}`[[1]]]))]
      traceTable[
        list(jj),
        r_i:=list(list(r_i[[1]] + traceTable[list(ii), `r_{ji}`[[1]]]))]
    }
  }
  setkey(traceTable, nodeId)
  traceTable[, nodeId:=NULL]

  traceTable[metaI$N+1]$omega_i[[1]] <- list("ND")
  traceTable[metaI$N+1]$Phi_i[[1]] <- list("ND")
  traceTable[metaI$N+1]$V_i[[1]] <- list("ND")
  traceTable[metaI$N+1]$V_1_i[[1]] <- list("ND")

  traceTable[metaI$N+1]$A_i[[1]] <- list("ND")
  traceTable[metaI$N+1]$b_i[[1]] <- list("ND")
  traceTable[metaI$N+1]$C_i[[1]] <- list("ND")
  traceTable[metaI$N+1]$d_i[[1]] <- list("ND")
  traceTable[metaI$N+1]$E_i[[1]] <- list("ND")
  traceTable[, ff_i:=as.list(f_i)]
  traceTable[, f_i:=NULL]
  setnames(traceTable, "ff_i", "f_i")
  traceTable[metaI$N+1]$f_i[[1]] <- list("ND")

  traceTable[metaI$N+1]$`L_{ji}`[[1]] <- list("ND")
  traceTable[metaI$N+1]$`m_{ji}`[[1]] <- list("ND")
  traceTable[, `rr_{ji}`:=as.list(`r_{ji}`)]
  traceTable[, `r_{ji}`:=NULL]
  setnames(traceTable, "rr_{ji}", "r_{ji}")
  traceTable[metaI$N+1]$`r_{ji}`[[1]] <- list("ND")


  traceTable[, `\\hat{X}_i`:=lapply(seq_len(.N), function(nodeId) {
    if(nodeId <= metaI$N) {
      "ND"
    # } else if(nodeId == metaI$N+1) {
    #   if(is.null(model$X0) || isTRUE(all(is.na(model$X0)))) {
    #     # set the root value to the one that maximizes the likelihood
    #     X0 <- try(solve(
    #       a=L_i[[nodeId]] + t(L_i[[nodeId]]),
    #       b = -m_i[[nodeId]]), silent = TRUE)
    #     if(inherits(X0, "try-error")) {
    #       X0 <- rep(NA_real_, PCMNumTraits(model))
    #     }
    #     X0
    #   } else {
    #     model$X0
    #   }
    } else {
      x <- rep(NaN, PCMNumTraits(model))
      if(any(is.finite(k_i[[nodeId]]))) {
        xHat <- try(
          solve(
            a = L_i[[nodeId]][k_i[[nodeId]], k_i[[nodeId]], drop = FALSE] +
              t(L_i[[nodeId]][k_i[[nodeId]],k_i[[nodeId]], drop = FALSE]),
            b = -m_i[[nodeId]][k_i[[nodeId]]]), silent = TRUE)
      }
      if(inherits(xHat, "try-error")) {
        x[k_i[[nodeId]]] <- rep(NA_real_, PCMNumTraits(model))
      } else {
        x[k_i[[nodeId]]] <- xHat
      }
      x
    }
  })]

  traceTable[, `\\ell\\ell_i`:=lapply(seq_len(.N), function(nodeId) {

    if(nodeId <= metaI$N) {
      # tip node
      "ND"
    } else {
      X0 <- `\\hat{X}_i`[[nodeId]]

      # if(nodeId == metaI$N+1) {
      #   # root node
      #   ll <- try(X0 %*% L_i[[nodeId]] %*% X0 + m_i[[nodeId]] %*% X0 + r_i[[nodeId]], silent = TRUE)
      #   if(inherits(ll, "try-error")) {
      #     ll <- NA_real_
      #   }
      #   ll
      # } else {
        # internal node
        ll <- try(
          X0[k_i[[nodeId]]] %*%
            L_i[[nodeId]][k_i[[nodeId]],k_i[[nodeId]], drop=FALSE] %*%
            X0[k_i[[nodeId]]] +

          m_i[[nodeId]][k_i[[nodeId]]] %*% X0[k_i[[nodeId]]] +

          r_i[[nodeId]],
          silent = TRUE)
        if(inherits(ll, "try-error")) {
          ll <- NA_real_
        }
        ll
      #}
    }
  })]

  traceTable
}

#' Quadratic polynomial parameters A, b, C, d, E, f for each node
#' @description An S3 generic function that has to be implemented for every
#'  model class. This function is called by \code{\link{PCMLik}}.
#' @inheritParams PCMLik
#' @export
PCMAbCdEf <- function(
  tree, model,
  SE = matrix(0.0, PCMNumTraits(model), PCMTreeNumTips(tree)),
  metaI=PCMInfo(NULL, tree, model, verbose = verbose),
  verbose = FALSE) {

  UseMethod("PCMAbCdEf", model)
}

#' @export
PCMAbCdEf.default <- function(
  tree, model,
  SE = matrix(0.0, PCMNumTraits(model), PCMTreeNumTips(tree)),
  metaI = PCMInfo(NULL, tree, model, verbose = verbose),
  verbose = FALSE) {

  threshold_SV <- metaI$PCMBase.Threshold.SV
  skip_singular <- metaI$PCMBase.Skip.Singular
  threshold_skip_singular <- metaI$PCMBase.Threshold.Skip.Singular
  threshold_eigval <- metaI$PCMBase.Threshold.EV

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

  omega <- array(NA_real_, dim=c(k, M))
  Phi <- array(NA_real_, dim=c(k, k, M))
  V <- array(NA_real_, dim=c(k, k, M))
  V_1 <- array(NA_real_, dim=c(k, k, M))


  # returned general form parameters
  A <- array(NA_real_, dim=c(k, k, M))
  b <- array(NA_real_, dim=c(k, M))
  C <- array(NA_real_, dim=c(k, k, M))
  d <- array(NA_real_, dim=c(k, M))
  E <- array(NA_real_, dim=c(k, k, M))
  f <- array(NA_real_, dim=c(M))

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

    # standard error
    if(i <= N) {
      VE <- metaI$VE[, , i]
    } else {
      VE <- matrix(0.0, metaI$k, metaI$k)
    }

    # length of edge leading to i
    ti <- tree$edge.length[edgeIndex]

    # present coordinates in parent and daughte nodes
    kj <- pc[,j]
    ki <- pc[,i]

    omega[,i] <- cond[[r[edgeIndex]]]$omega(ti, edgeIndex, metaI)
    Phi[,,i] <- cond[[r[edgeIndex]]]$Phi(ti, edgeIndex, metaI)
    V[,,i] <- cond[[r[edgeIndex]]]$V(ti, edgeIndex, metaI) + VE

    # Ensure that V[,,i] is symmetric. This is to prevent numerical errors
    # but can potentially hide a logical error in the model classes or the
    # parameters.
    V[,,i] <- 0.5*(V[,,i]+t(V[,,i]))

    # check that V[ki,ki,] is non-singular
    svdV = svd(matrix(V[ki,ki,i], sum(ki)), 0, 0)$d
    eigval <- eigen(V[ki,ki,i], symmetric = TRUE, only.values = TRUE)$values

    if(is.na(min(svdV)/max(svdV)) ||
       min(svdV)/max(svdV) < threshold_SV ||
       eigval[sum(ki)] < threshold_eigval) {
      singular[edgeIndex] <- TRUE
      if(!skip_singular || i <= N || ti > threshold_skip_singular ) {
        err <- paste0(
          "GaussianPCM.R:PCMAbCdEf.default:",i,":",
          " The matrix V for node ", i, " (branch length=", ti, ")",
          " is nearly singular: min(svdV)/max(svdV)=", min(svdV)/max(svdV),
          "; smallest eigenvalue = ", eigval[sum(ki)])
        stop(err)
      }
    }

    if(!singular[edgeIndex]) {
      # check V is positive definite
      if(sum(ki) == 1) {
        if(V[ki,ki,i] < 0) {
          err <- paste0(
            "GaussianPCM.R:PCMAbCdEf.default:",i,":",
            " The matrix V[ki,ki,i] for node ", i, " (sum(ki)=", sum(ki), ")",
            " is not positive definite, V[ki,ki,i]. Check the model parameters; V[ki,ki,i]=", V[ki,ki,i], ".")
          stop(err)
        }
      } else {
        if(!isSymmetric(V[ki,ki,i], tol = metaI$PCMBase.Tolerance.Symmetric)) {
          err <- paste0(
            "GaussianPCM.R:PCMAbCdEf.default:",i,":",
            " The matrix V[ki,ki,i] for node ", i, " (sum(ki)=", sum(ki), ")",
            " is not symmetric. It must be symmetric positive definite.")
          stop(err)
        }

        if(any(eigval<=0)) {
          err <- paste0(
            "GaussianPCM.R:PCMAbCdEf.default:",i,":",
            " The matrix V[ki,ki,i] for node ", i, " (sum(ki)=", sum(ki), ")",
            " is not positive definite. It must be symmetric positive definite; eigval=", toString(eigval), ".")
          stop(err)
        }
      }

      V_1[ki,ki,i] <- solve(as.matrix(V[ki,ki,i]))
      # TODO: V_1.slice(i)(ki, ki) = real(eigvec * diagmat(1/eigval) * eigvec.t());

      `%op%` <- if(sum(ki) > 1) `%*%` else `*`

      A[ki,ki,i] <- (-0.5*V_1[ki,ki,i])
      E[kj,ki,i] <- t(Phi[ki,kj,i]) %op% V_1[ki,ki,i]
      b[ki,i] <- V_1[ki,ki,i] %*% omega[ki,i]
      C[kj,kj,i] <- -0.5 * matrix(E[kj,ki,i], sum(kj), sum(ki)) %*% matrix(Phi[ki,kj,i], sum(ki), sum(kj))
      d[kj,i] <- -E[kj,ki,i] %op% omega[ki,i]
      f[i] <- -0.5 * (t(omega[ki,i]) %*% V_1[ki,ki,i] %*% omega[ki,i] + sum(ki)*log(2*pi) + log(det(as.matrix(V[ki,ki,i]))))
    }
  }

  list(omega=omega, Phi=Phi, V=V, V_1=V_1,
       A=A, b=b, C=C, d=d, E=E, f=f,
       singular = singular)
}

#' Quadratic polynomial parameters L, m, r
#'
#' @inheritParams PCMLik
#' @param root.only logical indicating whether to return the calculated values of L,m,r
#'  only for the root or for all nodes in the tree.
#' @return A list with the members A,b,C,d,E,f,L,m,r for all nodes in the tree or
#'   only for the root if root.only=TRUE.
#' @export
PCMLmr <- function(
  X, tree, model,
  SE = matrix(0.0, PCMNumTraits(model), PCMTreeNumTips(tree)),
  metaI = PCMInfo(
    X = X, tree = tree, model = model, SE = SE, verbose = verbose),
  root.only = TRUE, verbose = FALSE) {
  UseMethod("PCMLmr", metaI)
}

#' @export
PCMLmr.default <- function(
  X, tree, model,
  SE = matrix(0.0, PCMNumTraits(model), PCMTreeNumTips(tree)),
  metaI = PCMInfo(
    X = X, tree = tree, model = model, SE = SE, verbose = verbose),
  root.only = FALSE,
  verbose = FALSE
) {
  N <- metaI$N; M <- metaI$M; k <- metaI$k;

  edge <- tree$edge
  pc <- metaI$pc

  AbCdEf <- PCMAbCdEf(
    tree = tree, model = model, SE = SE, metaI = metaI,
    verbose = verbose)

  L <- array(NA_real_, dim=c(k, k, M))
  m <- array(NA_real_, dim=c(k, M))
  r <- array(NA_real_, dim=c(M))

  for(i in seq(N+1, M)) {
    L[pc[,i],pc[,i],i] <- 0.0
    m[pc[,i],i] <- 0.0
    r[i] <- 0.0
  }

  postorder <- rev(metaI$preorder)

  # avoid redundant calculation
  log2pi <- log(2*pi)

  for(edgeIndex in postorder) {
    # parent and daughter nodes
    j <- edge[edgeIndex, 1]
    i <- edge[edgeIndex, 2]
    # present coordinates
    kj <- pc[, j]
    ki <- pc[, i]

    if(i <= N) {
      # all es pointing to tips
      if(!AbCdEf$singular[edgeIndex]) {
        # ensure symmetry of L[,,i]
        L[,,i] <- 0.5 * (AbCdEf$C[,,i] + t(AbCdEf$C[,,i]))

        r[i] <- with(AbCdEf, t(X[ki,i]) %*% A[ki,ki,i] %*% X[ki,i] +
                       t(X[ki,i]) %*% b[ki,i] + f[i])

        m[kj,i] <- with(AbCdEf, d[kj,i] + matrix(E[kj,ki,i], sum(kj), sum(ki)) %*% X[ki,i])
      }
    } else {
      # edge pointing to internal nodes, for which all children
      # nodes have been visited
      if(!AbCdEf$singular[edgeIndex]) {
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

    # add up to parent
    L[,, j] <- L[,,j] + L[,,i]
    m[, j] <- m[,j] + m[,i]
    r[j] <- r[j] + r[i]
  }

  if(root.only) {
    list(L = L[,,N+1],
         m = m[,N+1],
         r = r[N+1])
  } else {
    c(AbCdEf[c("omega", "Phi", "V", "V_1", "A", "b", "C", "d", "E", "f")],
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
#' @param Sigma a numerical k x k matrix - neutral drift unit-time variance-covariance matrix.
#' @param Sigmae a numerical k x k matrix - environmental variance-covariance matrix.
#' @param Sigmaj is the variance matrix of the normal jump distribution (default is NULL).
#' @param xi a vector of 0's and 1's corresponding to each branch in the tree. A value of 1
#' indicates that a jump takes place at the beginning of the branch. This arugment is only
#' used if Sigmaj is not NULL. Default is NULL.
#' @param e_Ht a numerical k x k matrix - the result of the matrix exponential expm(-t*H).
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
  H, Sigma, Sigmae = NULL, Sigmaj = NULL, xi = NULL,
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
    if(!is.null(Sigmae) && metaI$edge[edgeIndex,2] <= metaI$N) {
      res <- res + Sigmae
    }
    Re(res)
  }
}

#' @importFrom mvtnorm rmvnorm
PCMCondRandom <- function(PCMCondObject, n=1, x0, t, edgeIndex, metaI, VE) {

  Mu <- PCMCondObject$omega(t, edgeIndex, metaI) +
    PCMCondObject$Phi(t, edgeIndex, metaI)%*%x0

  Sigma <- PCMCondObject$V(t, edgeIndex, metaI) + VE
  # ensure symmetry for Sigma:
  Sigma <- 0.5 * (Sigma + t(Sigma))

  rmvnorm(n = n,
          mean = Mu,
          sigma = Sigma)

}

#' @importFrom mvtnorm dmvnorm
PCMCondDensity <- function(PCMCondObject, x, x0, t, edgeIndex, metaI, VE, log=FALSE) {
  dmvnorm(x,
          mean = PCMCondObject$omega(t, edgeIndex, metaI) + PCMCondObject$Phi(t, edgeIndex, metaI)%*%x0,
          sigma = PCMCondObject$V(t, edgeIndex, metaI) + VE,
          log=log)
}
