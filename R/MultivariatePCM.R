#' @name PCMBase
#'
#' @title A General Framework for Gaussian Phylogenetic Comparative Models
#'
#' @description
#'
#'
NULL

#' Breadth-first tree traversal
#' @param tree a phylo object with possible singleton nodes (i.e. internal nodes with
#' one daughter node)
#' @return a vector of indices of edges in tree$edge in breadth-first order.
orderBreadthFirst <- function(tree) {
  # number of tips
  N <- length(tree$tip.label)

  # total number of nodes in the tree is the number of edges + 1 for the root
  M <- dim(tree$edge)[1]+1

  ordFrom <- order(tree$edge[,1])

  # we need the ordered edges in order to easily traverse all edges starting from
  # a given node
  iFrom <- match(1:M, tree$edge[ordFrom, 1])

  # the result is a vector of edge indices in the breadth-first search order
  res <- vector(mode='integer', length=M-1)

  # node-indices at the current level (start from the root)
  cn <- N+1
  j <- 1
  while(length(cn)>0) {
    cnNext <- c()
    for(n in cn) {
      # if not a tip
      if(n>N) {
        # indices in ordFrom of edges starting from the current node
        if(n<M) {
          es <- iFrom[n]:(iFrom[n+1]-1)
        } else {
          es <- iFrom[n]:(M-1)
        }
        jNext <- j+length(es)
        res[j:(jNext-1)] <- ordFrom[es]
        j <- jNext
        cnNext <- c(cnNext, tree$edge[ordFrom[es], 2])
      }
    }
    cn <- cnNext
  }
  res
}

#' Extract information for pruning a tree used as cache during likelihood
#' calculation
#'
#' @param tree a phylo object
#'
#' @details The function is very slow on strongly unbalanced trees due to the
#' slow vector append() operation in R.
#'
#' @return a list of objects used by the R implementation of mvlik().
#' @export
pruneTree <- function(tree) {
  # number of tips
  N <- length(tree$tip.label)
  # number of all nodes
  M <- nrow(tree$edge)+1

  # order the edge-indices in increasing index of ending node
  endingAt <- order(rbind(tree$edge, c(0, N+1))[, 2])

  edge <- tree$edge

  # count the number of non-visited children for each internal node
  nvc <- rep(0, M)

  # indices of parent node for edges that haven't still been gone through
  # initially, this are all edges
  ee1 <- edge[, 1]

  while(length(ee1)) {
    # For every element of (N+1):M its index in ee1 or NA
    matchp <- match((N+1):M, ee1)
    matchp <- matchp[!is.na(matchp)]
    # add one unvisited chlidren to each parent node's nvc
    nvc[ee1[matchp]] <- nvc[ee1[matchp]] + 1
    # remove the edges we've just gone through
    ee1 <- ee1[-matchp]
  }

  # start from the edges leading to tips
  nodesVector <- c()
  nodesIndex <- c(0)

  unVector <- c()
  # pointers to last unique indices (un) in unVector
  unIndex <- c(0)

  # internal or tip- nodes to which we are currently pointing, i.e. we are at
  # their parent-nodes and we are about to process the brances leading to them.
  nodes <- 1:N

  while(nodes[1] != N+1) {
    nodesIndex <- c(nodesIndex, nodesIndex[length(nodesIndex)]+length(nodes))
    nodesVector <- c(nodesVector, nodes)

    # indices of edges that end at nodes
    es <- endingAt[nodes]
    nodes <- c()

    while(length(es)>0) {
      # unique index of every edge ending at some of the nodes
      un <- match(unique(edge[es, 1]), edge[es, 1])
      # add these indices to unVector
      unVector <- c(unVector, un)
      # index of the last element of the current un in unVector
      unIndex <- c(unIndex, unIndex[length(unIndex)]+length(un))
      nvc[edge[es[un], 1]] <- nvc[edge[es[un], 1]] - 1
      nodes <- c(nodes, edge[es[un][nvc[edge[es[un], 1]] == 0], 1])
      es <- es[-un]
    }
  }
  list(# all raws from edge, times t and regimes must be accessed using indices
       # from the edingAt vector.
       endingAt=endingAt,

       nodesVector=nodesVector,
       nodesIndex=nodesIndex,
       nLevels=length(nodesIndex)-1,
       unVector=unVector,
       unIndex=unIndex)
}


#' For every node (root, internal or tip) in tree, build a logical vector of
#' length k with TRUE values for every present coordinate.
#' @param X numeric k x N matrix of observed values, with possible NA entries. The
#' columns in X are in the order of tree$tip.label
#' @param tree a phylo object
#' @param pruneI a list returned by the pruneTree function. Either leave this
#' as default or pass a previously computed pruneI for the same tree.
#' @return a k x M logical matrix
#' @export
presentCoordinates <- function(X, tree, pruneI=pruneTree(tree)) {
  edge <- tree$edge
  endingAt <- pruneI$endingAt
  nodesVector <- pruneI$nodesVector
  nodesIndex <- pruneI$nodesIndex
  nLevels <- pruneI$nLevels
  unVector <- pruneI$unVector
  unIndex <- pruneI$unIndex
  unJ <- 1

  N <- length(tree$tip.label)
  M <- nrow(edge)+1
  k <- dim(X)[1]


  if(getOption("PCMBase.Internal.PC.Full", TRUE)) {
    pc <- rep(TRUE, k*M)
    dim(pc) <- c(k, M)
    pc[, 1:N] <- !is.na(X[, 1:N])
  } else {
    pc <- rep(FALSE, k*M)
    dim(pc) <- c(k, M)

    for(i in 1:nLevels) {
      nodes <- nodesVector[(nodesIndex[i]+1):nodesIndex[i+1]]
      es <- endingAt[nodes]

      if(nodes[1] <= N) {
        # all es pointing to tips
        pc[, nodes] <- !is.na(X[, nodes])
      } else {
        # edges pointing to internal nodes, for which all children nodes have been
        # visited
        # here we do nothing
      }

      #update parent pifs
      while(length(es)>0) {
        un <- unVector[(unIndex[unJ]+1):unIndex[unJ+1]]
        unJ <- unJ+1
        pc[, edge[es[un], 1]] <- pc[, edge[es[un], 1]] | pc[, edge[es[un], 2]]
        es <- es[-un]
      }
    }
    if(any(rowSums(pc) == 0)) {
      stop("ERR:02001:PCMBase:MultivariatePCM.R:presentCoordinates:: Some tips have 0 present coordinates. Consider removing these tips.")
    }
  }
  pc
}

#' Generate trait data on a tree according to a multivariate stochastic
#' model with one or several regimes
#'
#' @param x0 a numeric vector of length k (the number of traits) specifying the #' trait values at the root of the tree.
#' @param tree a phylo object specifying a rooted tree.
#' @param model an S3 object specifying the model.
#' @param metaI a named list containg meta-information about the data and the
#' model.
#' A: numeric k x k matrix specifying the selection strengths of the OU process;
#' Theta numeric k-vector specifying the long-term mean of the OU process
#' @param Sigma numeric k x k matrix specifying the stochastic time-unit variance of the
#' OU process;
#'
#' @param A,Theta,Sigma parameters of the OU process:
#' A: a R x k x k array, where R is the number of regimes of the
#' OU process, k is the number of variables (traits), each A[r,,]
#' containing the matrix A for regime r;
#' Theta: a R x k matrix, row Theta[r, ] containing the long-term
#' mean theta for regime r;
#' Sigma: a R x k x k array, each Sigma[r,,] containing the
#' matrix Sigma for regime r;
#' Sigmae: a R x k matrix, each Sigmae[r, ] containing the environmental
#' variances for the k traits in regime r
#'
#' @param X0 numeric vector of length k indicating the root value
#'
#' @return a numeric M x k matrix of values at all nodes of the tree (root,
#' internal and tip), where M is the number such nodes: M=dim(tree$edge)[1]+1,
#' with indices from 1 to N=length(tree$tip.label) corresponding to tips, N+1
#' corresponding to the root and bigger than N+1 corresponding to internal nodes.
#'
#' @importFrom mvtnorm rmvnorm
#' @export
mvsim <- function(
  tree, model, X0,
  metaI = validateModel(tree = tree, model = model,
                        verbose = verbose),
  verbose = FALSE) {

  if(length(X0)!=metaI$k) {
    stop(paste('ERR:02002:PCMBase:MultivariatePCM.R:mvsim:: X0 must be of length', metaI$k, '.'))
  }

  values <- errors <- matrix(0, nrow=metaI$k, ncol=dim(tree$edge)[1]+1)
  values[, metaI$N + 1] <- X0

  ordBF <- orderBreadthFirst(tree)

  # create a list of random generator functions for each regime
  funMVCond <- lapply(1:metaI$R, function(r) {
    mvcond(tree, model = model, r = r, verbose = verbose)$mvr
  })

  for(e in ordBF) {
    values[, tree$edge[e, 2]] <-
      funMVCond[[metaI$regimes[e]]](
        n=1, x0 = values[, tree$edge[e,1]], t = tree$edge.length[e], e = e)
    if(!is.null(model$Sigmae)) {
      errors[, tree$edge[e, 2]] <-
        rmvnorm(1, rep(0, metaI$k),
                as.matrix(model$Sigmae[,, metaI$regimes[e]]))
    }
  }

  list(values=values, errors=errors)
}


# The specifics of every model are programmed in specifications of a
# few S3 generic functions:

#' An S3 generic function validating a model for a given tree.
#' @description This function should be implemented for each model.
#' @export
validateModel <- function(tree, model, verbose = FALSE) {
  UseMethod("validateModel", model)
}

#' Specify a phylogenetic comparative model for a given tree
#'
#' @param tree a phylo object
#' @param modelName a character specifying the name of the model. The
#' name of the model should be used as an S3 class of all model
#' instances.
#' @param k integer, specifying the number of traits;
#' @param R integer, specifying the number of regimes;
#' @param paramNames a list of characters specifying the names of
#'  the model parameters. The length of the list is used to
#'  determine the number of parameters.
#' @param paramStorageModes a list of characters of the same length
#'  as `paramNames` specifying the storage.mode (native type) for
#'  each parameter. Each element of this list should be one of
#'  "double", "integer", "logical", "character", "complex", "raw".
#' @param paramDims a list of integer vectors, or character
#'   strings evaluating as integer vectors specifying the
#'   dimension of each model parameter for all model regimes.
#'   This list should be the same length as `paramNames`. If an element of the
#'   list is a character string it must evaluate to a vector of positive
#'   integers, based on the objects accessible from the model parameters
#'   (including k, R and ...) and
#'   the symbols N and M, where N stays for the number of tips in the
#'   tree, and M denotes the number of nodes. For example, if
#'   there are 5 traits and 2 regimes, the dimension of the parameter Alpha of
#'   an OU model should be specified as c(5, 5, 2) or "c(5, 5, 2)" or
#'   "c(k, k, R)". If a parameter is a flat vector or an atom, for which the
#'   dimension is NULL, specify an integer corresponding to its length, e.g.
#'   the parameter specifying the presence of jumps in an OU model with jumps,
#'   should have a dimension "M-1" corresponding to the number of edges in the
#'   tree.
#' @param ... additional arguments used in the evaluation of
#'   `paramDims`.
#' @seealso \code{\link{storage.mode}}
#' @return an S3 object of class ModelSpecification.
#' @export
specifyModel <- function(tree, modelName,
                         k, R, regimesUnique,
                         paramNames,
                         paramStorageModes,
                         paramDims,
                         ...) {

  if(class(tree) != "phylo") {
    stop("ERR:02003:PCMBase:MultivariatePCM.R:specifyModel:: tree should be an object of S3 class 'phylo'.")
  }

  if(!is.character(modelName) & length(modelName) == 1) {
    stop("ERR:02004:PCMBase:MultivariatePCM.R:specifyModel:: modelName should be a character string.")
  }

  spec <- list(
    modelName = as.character(modelName)[1],
    k = as.integer(k),
    R = as.integer(R),
    regimesUnique = regimesUnique,
    N = as.integer(length(tree$tip.label)),
    M = as.integer(nrow(tree$edge) + 1),
    paramNames = lapply(paramNames, as.character),
    paramStorageModes = lapply(paramStorageModes, as.character)
  )

  if(! (storage.mode(regimesUnique) %in% c("integer", "character")) ) {
    stop("ERR:02005:PCMBase:MultivariatePCM.R:specifyModel:: the storage mode of regimesUnique should be either integer or character")
  }

  if( storage.mode(regimesUnique) == "integer" ) {
    if(!identical(regimesUnique, 1:max(regimesUnique)))
      stop("ERR:02006:PCMBase:MultivariatePCM.R:specifyModel:: when regimesUnique is integer it should be an integer vector from 1 to max(regimesUnique), so it can be used as index in an array.")
  }

  paramDims = with(
    c(spec, list(...)),
    lapply(paramDims, function(dim) {
      if(is.character(dim)) {
        as.integer(eval(parse(text = dim)))
      } else {
        as.integer(dim)
      }
    }))

  spec[["paramDims"]] <- paramDims

  if(spec[["k"]] <= 0) {
    stop(paste0("ERR:02007:PCMBase:MultivariatePCM.R:specifyModel:: k should be a positive integer but was ", spec[["k"]], "."))
  }
  if(spec[["R"]] <= 0) {
    stop(paste0("ERR:02008:PCMBase:MultivariatePCM.R:specifyModel:: R should be a positive integer but was", spec[["R"]], "."))
  }

  if(spec[["N"]] <= 0) {
    stop(paste0("ERR:02009:PCMBase:MultivariatePCM.R:specifyModel:: The tree should have at least one tip but N was ", spec[["N"]], "."))
  }

  if(spec[["M"]] <= 1) {
    stop(paste0("ERR:02010:PCMBase:MultivariatePCM.R:specifyModel:: The tree should have at least two nodes, that is, a minimum of one tip and one root-node but M was ", spec[["M"]], "."))
  }

  lapply(spec[["paramDims"]], function(dim) {
   if(any(is.na(dim)) | any(dim <= 0)) {
     stop(paste0("ERR:02011:PCMBase:MultivariatePCM.R:specifyModel:: paramDims should be a list of  positive integer vectors or a vector of characters, each of which can be parsed into a positive  integer vector, but was ", toString(paramDims)))
   }
  })

  nParams <- length(paramNames)

  if(length(paramStorageModes) != nParams) {
    stop(paste0("ERR:02012:PCMBase:MultivariatePCM.R:specifyModel:: paramStorageModes should be the same length as paramNames but its length was ", length(paramStorageModes), "."))
  }

  if(!all(paramStorageModes %in% c("character", "integer", "double", "complex", "logical"))) {
    stop('ERR:02013:PCMBase:MultivariatePCM.R:specifyModel:: all paramStorageModes should be in {"character", "integer", "double", "complex", "logical"}.')
  }

  if(length(paramDims) != nParams) {
    stop("ERR:02014:PCMBase:MultivariatePCM.R:specifyModel:: paramDims should be the same length as paramNames but its length was ", length(paramDims), ".")
  }

  class(spec) <- "ModelSpecification"
  spec
}

#' General function for validating a model
#' @param tree an object of S3 class "phylo"
#' @param model a model object of S3 class the name of the model
#' @param modelSpec an object of S3 class ModelSpecification
#' @param verbose a logical indicating if the log messges should be
#'   printed on the screen during validation.
#'
validateModelGeneral <- function(tree, model, modelSpec, verbose) {
  if(class(model)[[1]] != modelSpec$modelName) {
    stop(paste0("ERR:02020:PCMBase:MultivariatePCM.R:validateModelGeneral:: The model object should be of S3 class: ", modelSpec$modelName, " but was of class ", class(model)[[1]], "."))
  }
  nParams <- length(modelSpec$paramNames)

  if(!all(modelSpec$paramNames %in% names(model))) {
    stop(paste0("ERR:02021:PCMBase:MultivariatePCM.R:validateModelGeneral:: model should be a named list with elements ", toString(modelSpec$paramNames)))
  }

  # number of tips
  N <- length(tree$tip.label)

  # number of nodes including the root
  M <- nrow(tree$edge)+1
  if(verbose) {
    cat('M=', M, '\n')
  }

  if(N != modelSpec$N | M != modelSpec$M) {
    stop("ERR:02022:PCMBase:MultivariatePCM.R:validateModelGeneral:: the numbers of tips and nodes in tree do not match modelSpec$N and modelSpec$M.")
  }

  # number of regimes
  if(is.null(tree$edge.regime)) {
    if(!is.null(names(tree$edge.length))) {
      tree$edge.regime <- names(tree$edge.length)
      regimes <- unique(tree$edge.regime)
      R <- length(regimes)
    } else {
      regimes <- 1
      R <- 1
    }

  } else {
    regimes <- unique(tree$edge.regime)
    R <- length(regimes)
  }

  if(verbose) {
    cat('R=',R,'\n')
  }

  if(verbose) {
    cat("Regimes found in tree$edge.regime: ")
    print(regimes)
  }

  if(R==1) {
    regimes <- rep(1, length(tree$edge.length))
  } else {
    regimes <- match(tree$edge.regime, modelSpec$regimesUnique)
    if(any(is.na(regimes))) {
      stop(paste0("ERR:02023:PCMBase:MultivariatePCM.R:validateModelGeneral::: ",
                  " Some of the regimes in tree$edge.regime not found in",
                  "modelSpec$regimesUnique.\n",
                  "tree$edge.regime=", toString(tree$edge.regime), "\n",
                  "modelSpec$regimesUnique=", toString(modelSpec$regimesUnique)))
    }
  }

  for(i in 1:nParams) {
    paramName_i <- modelSpec$paramNames[[i]]
    storageMode_i <- modelSpec$paramStorageModes[[i]]
    paramDim_i <- modelSpec$paramDims[[i]]
    if(storage.mode(model[[paramName_i]]) != storageMode_i) {
      stop(paste0(
        "ERR:02024:PCMBase:MultivariatePCM.R:validateModelGeneral:: ",
        paramName_i, " should have a storage mode ", storageMode_i,
        " but has storageMode", storage.mode(model[[paramName_i]]), "."))
    }
    if(length(paramDim_i) == 1) {
      # a vector or an atom
      if(length(model[[paramName_i]]) != paramDim_i) {
        stop(paste0(
          "ERR:02025:PCMBase:MultivariatePCM.R:validateModelGeneral:: ",
          paramName_i, " should be a vector of length ", paramDim_i,
          " but has length ", length(model[[paramName_i]])
        ))
      }
    } else {
      if(!identical(dim(model[[paramName_i]]), paramDim_i)) {
        stop(paste0(
          "ERR:02024:PCMBase:MultivariatePCM.R:validateModelGeneral:: ",
          paramName_i, " should have dimension ", toString(paramDim_i),
          " but has dimension ", dim(model[[paramName_i]]), "."
        ))
      }
    }
  }

  list(N = modelSpec$N, M = modelSpec$M, R = modelSpec$R, k=modelSpec$k,
       regimes = regimes, regimesUnique = modelSpec$regimesUnique)
}

#' Multivariate normal distribution for a given tree and model
#' @export
mvcond <- function(tree, model, metaI, r=1, verbose = FALSE) {
  UseMethod("mvcond", model)
}

#' Quadratic polynomial parameters A, b, C, d, E, f for each node
#' @description A generic function specified for every model class. This
#' function is called by mvlik.
#' @export
AbCdEf <- function(tree, model, metaI, pc, verbose = FALSE) {
  UseMethod("AbCdEf", model)
}

#' Quadratic polynomial parameters L, m, r
#' @export
Lmr <- function(X, tree, model,
                metaI = validateModel(tree, model, verbose = verbose),
                pruneI = pruneTree(tree),
                pc = presentCoordinates(X, tree, pruneI),
                root.only = TRUE, verbose = FALSE) {
  UseMethod("Lmr", pruneI)
}

#' Defaut R-implementation of Lmr
#' @details This funciton is not a generic S3 implementation of Lmr - it needs
#' additional raguments and is designed to be called by mvlik. It can still
#' be called by the end-user for debugging purpose.
#' @return A list with the members A,b,C,d,E,f,L,m,r for all nodes in the tree.
#' @export
Lmr.default <- function(X, tree, model, metaI=validateModel(tree, model),
                       pruneI = pruneTree(tree),
                       pc = presentCoordinates(X, tree, pruneI),
                       root.only = FALSE,
                       verbose = FALSE) {
  unJ <- 1

  N <- metaI$N; M <- metaI$M; k <- metaI$k;

  edge <- tree$edge
  endingAt <- pruneI$endingAt
  nodesVector <- pruneI$nodesVector
  nodesIndex <- pruneI$nodesIndex
  nLevels <- pruneI$nLevels
  unVector <- pruneI$unVector
  unIndex <- pruneI$unIndex

  threshold_SV <- getOption("PCMBase.Threshold.SV", 1e-6)

  L <- array(0, dim=c(k, k, M))
  m <- array(0, dim=c(k, M))
  r <- array(0, dim=c(M))


  AbCdEf <- AbCdEf(tree = tree, model = model,
                   metaI = metaI,
                   pc = pc, verbose = verbose)


  # avoid redundant calculation
  log2pi <- log(2*pi)

  for(level in 1:nLevels) {
    nodes <- nodesVector[(nodesIndex[level]+1):nodesIndex[level+1]]
    es <- endingAt[nodes]

    if(nodes[1] <= N) {
      # all es pointing to tips
      #L[nodes,,] <- AbCdEf$C[nodes,,]

      for(e in es) {
        # parent and daughter nodes
        j <- edge[e, 1]; i <- edge[e, 2];
        # present coordinates
        kj <- pc[, j]; ki <- pc[, i];


        # check that V[ki,ki,] is non-singular
        svdV = svd(matrix(AbCdEf$V[ki,ki,i], sum(ki)), 0, 0)$d
        if(min(svdV)/max(svdV) < threshold_SV) {
          err <- paste0(
            "ERR:02031:PCMBase:MultivariatePCM.R:Lmr.default:",i,":",
            " The matrix V for node ", i,
            " is nearly singular: min(svdV)/max(svdV)=", min(svdV)/max(svdV),
            ". Check the model parameters and the length of the branch",
            " leading to the node. For details on this error, read the User Guide.")
          stop(err)
        }

        # ensure symmetry of L[,,i]
        L[,,i] <- 0.5 * (AbCdEf$C[,,i] + t(AbCdEf$C[,,i]))

        r[i] <- with(AbCdEf, t(X[ki,i]) %*% A[ki,ki,i] %*% X[ki,i] +
                       t(X[ki,i]) %*% b[ki,i] + f[i])

        m[kj,i] <- with(AbCdEf, d[kj,i] + matrix(E[kj,ki,i], sum(kj), sum(ki)) %*% X[ki,i])

        #logDetV[i] <- with(AbCdEf, det(-2*(A[i,ki,ki]+L[i,ki,ki])))

        #K <- K + sum(ki)
      }
    } else {
      # edges pointing to internal nodes, for which all children
      # nodes have been visited
      for(e in es) {
        # parent and daughter nodes
        j <- edge[e, 1]; i <- edge[e, 2];
        # present coordinates
        kj <- pc[, j]; ki <- pc[, i];

        # check that V[i,ki,ki] is non-singular
        svdV = svd(matrix(AbCdEf$V[ki,ki,i], sum(ki)), 0, 0)$d
        if(min(svdV)/max(svdV) < threshold_SV) {
          err <- paste0(
            "ERR:02031:PCMBase:MultivariatePCM.R:Lmr.default:",i,":",
            " The matrix V for node ", i,
            " is nearly singular: min(svdV)/max(svdV)=", min(svdV)/max(svdV),
            ", det(V)=", det(matrix(AbCdEf$V[i,ki,ki], sum(ki))),
            ". Check the model parameters and the length of the branch",
            " leading to the node. For details on this error, read the User Guide.")
          stop(err)
        }



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


#' Multivariate likelihood calculation
#' @export
mvlik <- function(X, tree, model,
                metaI =validateModel(tree, model, verbose = verbose),
                pruneI = pruneTree(tree),
                pc = presentCoordinates(X, tree),
                log = TRUE,
                verbose = FALSE,
                debug = FALSE) {

  # will change this value if there is no error
  value.NA <- getOption("PCMBase.Value.NA", as.double(NA))

  Lmr <- try(Lmr(X, tree, model, metaI, pruneI, pc, verbose = verbose, root.only = TRUE),
             silent = TRUE)

  if(class(Lmr) == "try-error") {
    errL <- parseErrorMessage(Lmr)
    if(is.null(errL)) {
      err <- paste0("ERR:02041:PCMBase:MultivariatePCM.R:mvlik:: There was a problem calculating the coefficients L,m,r. Error message from call to Lmr: ", Lmr)
      errL <- parseErrorMessage(err)
    } else {
      err <- Lmr
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

  } else if(is.list(Lmr)) {

    L_root <- Lmr$L
    m_root <- Lmr$m
    r_root <- Lmr$r

    if(is.null(L_root) | is.null(m_root) | is.null(r_root)) {
      err <- paste0("ERR:02042:PCMBase:MultivariatePCM.R:mvlik:: The list returned by Lmr did not contain elements 'L', 'm' and 'r'.")
      if(getOption("PCMBase.Errors.As.Warnings", TRUE)) {
        warning(err)
      } else {
        stop(err)
      }

      errL <- parseErrorMessage(err)

      X0 <- model$X0
      attr(value.NA, 'X0') <- X0
      attr(value.NA, "error") <- errL

      return(value.NA)

    }

    if(is.null(model$X0)) {
      # set the root value to the one that maximizes the likelihood
      X0 <- try(solve(a=L_root + t(L_root), b = -m_root), silent = TRUE)
      if(class(X0) == "try-error") {
        err <- paste0(
          "ERR:02043:PCMBase:MultivariatePCM.R:mvlik:: There was a problem calculating X0 from the coefficients L,m,r. ", "L=", toString(L), "; m=", toString(m), "; r=", r,
          ". Error message from call to solve(a=L_root + t(L_root), b = -m_root):", X0)

        errL <- parseErrorMessage(err)
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
        "ERR:02044:PCMBase:MultivariatePCM.R:mvlik:: There was a problem calculating loglik from X0 and the coefficients L,m,r. ", "X0=", toString(X0), "L=", toString(L), "; m=", toString(m), "; r=", r,
        ". Error message from call to X0 %*% L_root %*% X0 + m_root %*% X0 + r_root:", loglik)

      errL <- parseErrorMessage(err)
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
        "ERR:02045:PCMBase:MultivariatePCM.R:mvlik:: There was a problem calculating value from loglik=", toString(loglik), ". Error message from call to as.vector(if(log) loglik else exp(loglik)):", value)

      errL <- parseErrorMessage(err)
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
        "ERR:02046:PCMBase:MultivariatePCM.R:mvlik:: There was a possible numerical problem, e.g. division of 0 by 0 when calculating the likelihood. value=", toString(value), "; calculated loglik=", toString(loglik), ". No error message was returned from the call to Lmr. Check for runtime warnings.")

      errL <- parseErrorMessage(err)
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

#' Extract error information from a formatted error message.
#' @param x character string representing the error message.
#' @description The function searches x for a pattern matching the format
#' 'ERR:5-digit-code:project-name:source-file:error-specifics:'. Specifically it
#' searches for a regular expression pattern "ERR:[0-9]+:[^:]+:[^:]+:[^:]+:[^:]*:".
#' @return a named list with the parsed error information or NULL, if no match
#' was found.
#' @export
parseErrorMessage <- function(x) {
  res <- try({
    if(is.character(x)) {
      code <- regmatches(x, regexpr("ERR:[0-9]+:[^:]+:[^:]+:[^:]+:[^:]*:", x))
      if(length(code) > 0) {
        code <- code[1]
        codeL <- strsplit(code, split=":")[[1]]
        list(
          type = codeL[1],
          icode = codeL[2],
          project = codeL[3],
          file = codeL[4],
          fun = codeL[5],
          info = codeL[6],
          code = code,
          msg = x
        )
      } else {
        NULL
      }
    } else {
      NULL
    }
  }, silent = TRUE)

  if(class(res)=="try-error") {
    NULL
  } else {
    res
  }
}
