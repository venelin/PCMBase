#' @name PCMBase
#'
#' @title A General Framework for Gaussian Phylogenetic Comparative Models
#'
#' @description
#'
#'
NULL

#' Simulation of a phylogenetic comparative model on a tree
#'
#' @description Generate trait data on a tree according to a multivariate stochastic
#' model with one or several regimes
#'
#' @param X0 a numeric vector of length k (the number of traits) specifying the
#' trait values at the root of the tree.
#' @param tree a phylo object specifying a rooted tree.
#' @param model an S3 object specifying the model (see Details).
#' @param metaI a named list containg meta-information about the data and the
#' model.
#'
#' @details Internally, this function uses the \code{\link{PCMCond}} iimplementation
#'  for the given model class.
#'
#' @return a list with two members as follows:
#' \item{values}{numeric M x k matrix of values at all nodes of the tree (root,
#' internal and tip), where M is the number of nodes: M=dim(tree$edge)[1]+1,
#' with indices from 1 to N=length(tree$tip.label) corresponding to tips, N+1
#' corresponding to the root and bigger than N+1 corresponding to internal nodes.}
#' \item{errors}{numeric M x k matrix of errors at all nodes of the tree (root,
#' internal and tip). Note that these errors are only simulated if the model contains
#' a k x k x R matrix member with the name 'Sigmae'. In this case the errors are
#' simulated as a random Gaussian noise with 0 mean and variance covariance matrix
#' defined by Sigmae. Note also that the errors on the root- and the internal nodes
#' are not used as input for the values/errors at their descending nodes, because
#' they are treated as non-heritable components or measurement error.}
#' The function will fail in case that the length of the argument vector X0 differs
#' from the number of traits specified in \code{metaI$k}. Error message: "ERR:02002:PCMBase:MultivariatePCM.R:PCMSim:: X0 must be of length ...".
#'
#' @importFrom mvtnorm rmvnorm
#' @seealso \code{\link{PCMLik}} \code{\link{PCMValidate}} \code{\link{PCMCond}}
#' @export
PCMSim <- function(
  tree, model, X0,
  metaI = PCMValidate(tree = tree, model = model,
                        verbose = verbose),
  verbose = FALSE) {

  if(length(X0)!=metaI$k) {
    stop(paste('ERR:02002:PCMBase:MultivariatePCM.R:PCMSim:: X0 must be of length', metaI$k, '.'))
  }

  values <- errors <- matrix(0, nrow=metaI$k, ncol=dim(tree$edge)[1]+1)
  values[, metaI$N + 1] <- X0

  ordBF <- PCMBreadthFirstOrder(tree)

  # create a list of random generator functions for each regime
  funMVCond <- lapply(1:metaI$R, function(r) {
    PCMCond(tree, model = model, r = r, verbose = verbose)$random
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

#' Likelihood of a multivariate Gaussian phylogenetic comparative model with non-interacting lineages
#'
#' @description The likelihood of a PCM represets the probability density function
#'   of observed trait values (data) at the tips of a tree given the tree and
#'   the model parameters. Seen as a function of the model parameters, the
#'   likelihood is used to fit the model to the observed trait data and the
#'   phylogenetic tree (which is typically inferred from another sort of data, such
#'   as an alignment of genetic sequences for the species at the tips of the tree).
#'   The \code{\link{PCMLik}} function
#'   provides a common interface for calculating the (log-)likelihood of different
#'   PCMs.
#'   Below we denote by N the number of tips, by M the total number of nodes in the
#'   tree including tips, internal and root node, and by k - the number of traits.
#'
#' @param X a \code{k x N} numerical matrix with possible \code{NA} entries. Each
#'   column of X contains the measured trait values for one species (tip in tree).
#'   Depending on the value of the global option "PCMBase.Internal.PC.Full" \code{NA}
#'   entries are interpreted either as missing measurements or non-existing traits
#'   (see \code{\link{PCMPresentCoordinates}}).
#' @param tree a phylo object with N tips.
#' @param model an S3 object specifying both, the model type (class, e.g. "OU") as
#'   well as the concrete model parameter values at which the likelihood is to be
#'   calculated (see also Details).
#' @param metaI a list returned from a call to \code{PCMValidate(tree, model)},
#'   containing meta-data such as N, M and k.
#' @param pruneI a named list containing cached preprocessing data for the tree used
#'   to perform post-order traversal (pruning). By default, this is created
#'   using \code{PCMPruningOrder(tree)}. This will use the default R-implementation of the
#'   likelihood calculation, which is based on the default R-implementation of the
#'   function \code{\link{PCMLmr}} (\code{\link{PCMLmr.default}}) and the S3 specification
#'   of the function
#'   \code{\link{PCMAbCdEf}} for the given model (a function called \code{PCMAbCdEf.MODEL},
#'   where model is the name of the mode and the class attribute of the \code{model}
#'   argument. For a different implementation of the function \code{\link{PCMLmr}},
#'   provide an S3 object for which an S3 specification has been written (see Details
#'   and example section).
#' @param pc a logical k x M matrix returned by \code{\link{PCMPresentCoordinates}}.
#' @param log logical indicating whether a log-liklehood should be calculated. Default
#'  is TRUE.
#' @param verbose logical indicating if some debug-messages should printed.
#'
#' @return a numerical value with named attributes as follows:
#' \item{X0}{A numerical vector of length k specifying the value at the root for which
#' the likelihood value was calculated. If the model contains a member called X0, this
#' vector is used; otherwise the value of X0 maximizing the likelihood for the given
#' model parameters is calculated by maximizing the quadratic polynomial
#' 'X0 * L_root * X0 + m_root * X0 + r_root'.}
#' \item{error}{A named list containing error information if a numerical or other
#' logical error occured during likelihood calculation (this is a list returned by
#'  \code{\link{PCMParseErrorMessage}}.}
#'  If an error occured during likelihood calculation, the default behavior is to
#'  return NA with a non-NULL error attribute. This behavior can be changed in
#'  using global options:
#'  \item{"PCMBase.Value.NA"}{Allows to specify a different NA value such as \code{-Inf} or \code{-1e20} which can be used in compbination with \code{log = TRUE} when
#'   using \code{optim} to maximize the log-likelihood;}
#'  \item{"PCMBase.Errors.As.Warnings"}{Setting this option to FALSE will cause any
#'  error to result in calling the \code{\link{stop}} R-base function. If not caught
#'  in a \code{\link{tryCatch}}, this will cause the inference procedure to abort at the occurence of a numerical error. By default, this option is set to TRUE, which
#'  means that \code{getOption("PCMBase.Value.NA", as.double(NA))} is returned with
#'  an error attribute and a warning is issued.}
#'
#' @details For efficiency, the arguments \code{metaI}, \code{pruneI} and \code{pc}
#'   can be provided explicitly, because these are not supposed to change during a
#'   model inference procedure such as likelihood maximization (see example).
#'
#' @seealso \code{\link{PCMValidate}} \code{\link{PCMPresentCoordinates}} \code{\link{PCMPruningOrder}} \code{\link{PCMAbCdEf}} \code{\link{PCMLmr}} \code{\link{PCMSim}} \code{\link{PCMCond}} \code{\link{PCMParseErrorMessage}}
#' @export
PCMLik <- function(X, tree, model,
                  metaI =PCMValidate(tree, model, verbose = verbose),
                  pruneI = PCMPruningOrder(tree),
                  pc = PCMPresentCoordinates(X, tree),
                  log = TRUE,
                  verbose = FALSE) {

  # will change this value if there is no error
  value.NA <- getOption("PCMBase.Value.NA", as.double(NA))

  PCMLmr <- try(PCMLmr(X, tree, model, metaI, pruneI, pc, verbose = verbose, root.only = TRUE),
             silent = TRUE)

  if(class(PCMLmr) == "try-error") {
    errL <- PCMParseErrorMessage(PCMLmr)
    if(is.null(errL)) {
      err <- paste0("ERR:02041:PCMBase:MultivariatePCM.R:PCMLik:: There was a problem calculating the coefficients L,m,r. Error message from call to PCMLmr: ", PCMLmr)
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
      err <- paste0("ERR:02042:PCMBase:MultivariatePCM.R:PCMLik:: The list returned by PCMLmr did not contain elements 'L', 'm' and 'r'.")
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

    if(is.null(model$X0)) {
      # set the root value to the one that maximizes the likelihood
      X0 <- try(solve(a=L_root + t(L_root), b = -m_root), silent = TRUE)
      if(class(X0) == "try-error") {
        err <- paste0(
          "ERR:02043:PCMBase:MultivariatePCM.R:PCMLik:: There was a problem calculating X0 from the coefficients L,m,r. ", "L=", toString(L), "; m=", toString(m), "; r=", r,
          ". Error message from call to solve(a=L_root + t(L_root), b = -m_root):", X0)

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
        "ERR:02044:PCMBase:MultivariatePCM.R:PCMLik:: There was a problem calculating loglik from X0 and the coefficients L,m,r. ", "X0=", toString(X0), "L=", toString(L), "; m=", toString(m), "; r=", r,
        ". Error message from call to X0 %*% L_root %*% X0 + m_root %*% X0 + r_root:", loglik)

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
        "ERR:02045:PCMBase:MultivariatePCM.R:PCMLik:: There was a problem calculating value from loglik=", toString(loglik), ". Error message from call to as.vector(if(log) loglik else exp(loglik)):", value)

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
        "ERR:02046:PCMBase:MultivariatePCM.R:PCMLik:: There was a possible numerical problem, e.g. division of 0 by 0 when calculating the likelihood. value=", toString(value), "; calculated loglik=", toString(loglik), ". No error message was returned from the call to PCMLmr. Check for runtime warnings.")

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


#' Breadth-first tree traversal
#' @param tree a phylo object with possible singleton nodes (i.e. internal nodes with
#' one daughter node)
#' @return a vector of indices of edges in tree$edge in breadth-first order.
PCMBreadthFirstOrder <- function(tree) {
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
#' @return a list of objects used by the R implementation of PCMLik().
#' @export
PCMPruningOrder <- function(tree) {
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


#' Determine which traits are present (active) on each node of the tree
#'
#' @description For every node (root, internal or tip) in tree, build a logical
#' vector of length k with TRUE values for every present coordinate. Non-present
#' coordinates arize from NA-values in the trait data. These can occur in two cases:
#' \describe{
#' \item{Missing measurements for some traits at some tips:}{the present coordinates are FALSE for the corresponding
#' tip and trait, but are full for all traits at all internal and root nodes.}
#' \item{non-existent traits for some species:}{the FALSE present coordinates propagate towards the parent
#' nodes - an internal or root node will have a present coordinate set to FALSE
#' for a given trait, if all of its descendants have this coordinate set to FALSE.}
#' }
#' These two cases have different effect on the likelihood calculation: missing
#' measurements are integrated out at the parent nodes; while non-existent traits
#' are treated as reduced dimensionality of the vector at the parent node.
#' The PCMBase package allows to specify the treatment of NA trait values using the
#' global option "PCMBase.Internal.PC.Full", which is set to TRUE by default. This
#' setting corresponds to assuming that all traits exist for all tips and NA values
#' are integrated out during likelihood calculation. Setting this option to FALSE
#' corresponds to assuming all NA values as non-existent.
#'
#' @param X numeric k x N matrix of observed values, with possible NA entries. The
#' columns in X are in the order of tree$tip.label
#' @param tree a phylo object
#' @param pruneI a list returned by the PCMPruningOrder function. Either leave this
#' as default or pass a previously computed pruneI for the same tree.
#' @return a k x M logical matrix which can be passed as a pc argument to the PCMLik
#' function. The function fails in case when all traits are NAs for some of the tips.
#' In that case an error message is issued "ERR:02001:PCMBase:MultivariatePCM.R:PCMPresentCoordinates:: Some tips have 0 present coordinates. Consider removing these tips.".
#' @seealso \code{\link{PCMLik}}
#' @export
PCMPresentCoordinates <- function(X, tree, pruneI=PCMPruningOrder(tree)) {
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
      stop("ERR:02001:PCMBase:MultivariatePCM.R:PCMPresentCoordinates:: Some tips have 0 present coordinates. Consider removing these tips.")
    }
  }
  pc
}


# The specifics of every model are programmed in specifications of a
# few S3 generic functions:

#' Validating a phylogenetic comparative model
#'
#' @description An S3 generic function validating a model for a given tree. This
#'   function should be implemented for each model class.
#' @param tree a phylo object
#' @param model an S3 object specifying a model (a list with a non-null class attribute).
#' @param verbose a logical specifying if some information should be printed on
#' the console.
#' @return a logical value
#' @export
PCMValidate <- function(tree, model, verbose = FALSE) {
  UseMethod("PCMValidate", model)
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
PCMSpecify <- function(tree, modelName,
                         k, R, regimesUnique,
                         paramNames,
                         paramStorageModes,
                         paramDims,
                         ...) {

  if(class(tree) != "phylo") {
    stop("ERR:02003:PCMBase:MultivariatePCM.R:PCMSpecify:: tree should be an object of S3 class 'phylo'.")
  }

  if(!is.character(modelName) & length(modelName) == 1) {
    stop("ERR:02004:PCMBase:MultivariatePCM.R:PCMSpecify:: modelName should be a character string.")
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
    stop("ERR:02005:PCMBase:MultivariatePCM.R:PCMSpecify:: the storage mode of regimesUnique should be either integer or character")
  }

  if( storage.mode(regimesUnique) == "integer" ) {
    if(!identical(regimesUnique, 1:max(regimesUnique)))
      stop("ERR:02006:PCMBase:MultivariatePCM.R:PCMSpecify:: when regimesUnique is integer it should be an integer vector from 1 to max(regimesUnique), so it can be used as index in an array.")
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
    stop(paste0("ERR:02007:PCMBase:MultivariatePCM.R:PCMSpecify:: k should be a positive integer but was ", spec[["k"]], "."))
  }
  if(spec[["R"]] <= 0) {
    stop(paste0("ERR:02008:PCMBase:MultivariatePCM.R:PCMSpecify:: R should be a positive integer but was", spec[["R"]], "."))
  }

  if(spec[["N"]] <= 0) {
    stop(paste0("ERR:02009:PCMBase:MultivariatePCM.R:PCMSpecify:: The tree should have at least one tip but N was ", spec[["N"]], "."))
  }

  if(spec[["M"]] <= 1) {
    stop(paste0("ERR:02010:PCMBase:MultivariatePCM.R:PCMSpecify:: The tree should have at least two nodes, that is, a minimum of one tip and one root-node but M was ", spec[["M"]], "."))
  }

  lapply(spec[["paramDims"]], function(dim) {
   if(any(is.na(dim)) | any(dim <= 0)) {
     stop(paste0("ERR:02011:PCMBase:MultivariatePCM.R:PCMSpecify:: paramDims should be a list of  positive integer vectors or a vector of characters, each of which can be parsed into a positive  integer vector, but was ", toString(paramDims)))
   }
  })

  nParams <- length(paramNames)

  if(length(paramStorageModes) != nParams) {
    stop(paste0("ERR:02012:PCMBase:MultivariatePCM.R:PCMSpecify:: paramStorageModes should be the same length as paramNames but its length was ", length(paramStorageModes), "."))
  }

  if(!all(paramStorageModes %in% c("character", "integer", "double", "complex", "logical"))) {
    stop('ERR:02013:PCMBase:MultivariatePCM.R:PCMSpecify:: all paramStorageModes should be in {"character", "integer", "double", "complex", "logical"}.')
  }

  if(length(paramDims) != nParams) {
    stop("ERR:02014:PCMBase:MultivariatePCM.R:PCMSpecify:: paramDims should be the same length as paramNames but its length was ", length(paramDims), ".")
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
#' @return a named list
PCMValidateGeneral <- function(tree, model, modelSpec, verbose) {
  if(class(model)[[1]] != modelSpec$modelName) {
    stop(paste0("ERR:02020:PCMBase:MultivariatePCM.R:PCMValidateGeneral:: The model object should be of S3 class: ", modelSpec$modelName, " but was of class ", class(model)[[1]], "."))
  }
  nParams <- length(modelSpec$paramNames)

  if(!all(modelSpec$paramNames %in% names(model))) {
    stop(paste0("ERR:02021:PCMBase:MultivariatePCM.R:PCMValidateGeneral:: model should be a named list with elements ", toString(modelSpec$paramNames)))
  }

  # number of tips
  N <- length(tree$tip.label)

  # number of nodes including the root
  M <- nrow(tree$edge)+1
  if(verbose) {
    cat('M=', M, '\n')
  }

  if(N != modelSpec$N | M != modelSpec$M) {
    stop("ERR:02022:PCMBase:MultivariatePCM.R:PCMValidateGeneral:: the numbers of tips and nodes in tree do not match modelSpec$N and modelSpec$M.")
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
      stop(paste0("ERR:02023:PCMBase:MultivariatePCM.R:PCMValidateGeneral::: ",
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
        "ERR:02024:PCMBase:MultivariatePCM.R:PCMValidateGeneral:: ",
        paramName_i, " should have a storage mode ", storageMode_i,
        " but has storageMode", storage.mode(model[[paramName_i]]), "."))
    }
    if(length(paramDim_i) == 1) {
      # a vector or an atom
      if(length(model[[paramName_i]]) != paramDim_i) {
        stop(paste0(
          "ERR:02025:PCMBase:MultivariatePCM.R:PCMValidateGeneral:: ",
          paramName_i, " should be a vector of length ", paramDim_i,
          " but has length ", length(model[[paramName_i]])
        ))
      }
    } else {
      if(!identical(dim(model[[paramName_i]]), paramDim_i)) {
        stop(paste0(
          "ERR:02024:PCMBase:MultivariatePCM.R:PCMValidateGeneral:: ",
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
#' @description An S3 generic function that has to be implemented for every
#'  model class.
#' @inheritParams PCMLik
#' @param r an integer specifying a model regime
#' @return a list with the following members:
#' \item{random}{}
#' \item{density}{}
#' @export
PCMCond <- function(tree, model, metaI, r=1, verbose = FALSE) {
  UseMethod("PCMCond", model)
}

#' Quadratic polynomial parameters A, b, C, d, E, f for each node
#' @description An S3 generic function that has to be implemented for every
#'  model class. This function is called by \code{\link{PCMLik}}.
#' @inheritParams PCMLik
#' @export
PCMAbCdEf <- function(tree, model, metaI, pc, verbose = FALSE) {
  UseMethod("PCMAbCdEf", model)
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
                metaI = PCMValidate(tree, model, verbose = verbose),
                pruneI = PCMPruningOrder(tree),
                pc = PCMPresentCoordinates(X, tree, pruneI),
                root.only = TRUE, verbose = FALSE) {
  UseMethod("PCMLmr", pruneI)
}

#' Defaut R-implementation of PCMLmr
#' @details This funciton is not a generic S3 implementation of PCMLmr - it needs
#' additional raguments and is designed to be called by PCMLik. It can still
#' be called by the end-user for debugging purpose.
#' @inheritParams PCMLik
#' @return A list with the members A,b,C,d,E,f,L,m,r for all nodes in the tree or
#'   only for the root if root.only=TRUE.
#' @export
PCMLmr.default <- function(X, tree, model, metaI=PCMValidate(tree, model),
                       pruneI = PCMPruningOrder(tree),
                       pc = PCMPresentCoordinates(X, tree, pruneI),
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


  PCMAbCdEf <- PCMAbCdEf(tree = tree, model = model,
                   metaI = metaI,
                   pc = pc, verbose = verbose)


  # avoid redundant calculation
  log2pi <- log(2*pi)

  for(level in 1:nLevels) {
    nodes <- nodesVector[(nodesIndex[level]+1):nodesIndex[level+1]]
    es <- endingAt[nodes]

    if(nodes[1] <= N) {
      # all es pointing to tips
      #L[nodes,,] <- PCMAbCdEf$C[nodes,,]

      for(e in es) {
        # parent and daughter nodes
        j <- edge[e, 1]; i <- edge[e, 2];
        # present coordinates
        kj <- pc[, j]; ki <- pc[, i];


        # check that V[ki,ki,] is non-singular
        svdV = svd(matrix(PCMAbCdEf$V[ki,ki,i], sum(ki)), 0, 0)$d
        if(min(svdV)/max(svdV) < threshold_SV) {
          err <- paste0(
            "ERR:02031:PCMBase:MultivariatePCM.R:PCMLmr.default:",i,":",
            " The matrix V for node ", i,
            " is nearly singular: min(svdV)/max(svdV)=", min(svdV)/max(svdV),
            ". Check the model parameters and the length of the branch",
            " leading to the node. For details on this error, read the User Guide.")
          stop(err)
        }

        # ensure symmetry of L[,,i]
        L[,,i] <- 0.5 * (PCMAbCdEf$C[,,i] + t(PCMAbCdEf$C[,,i]))

        r[i] <- with(PCMAbCdEf, t(X[ki,i]) %*% A[ki,ki,i] %*% X[ki,i] +
                       t(X[ki,i]) %*% b[ki,i] + f[i])

        m[kj,i] <- with(PCMAbCdEf, d[kj,i] + matrix(E[kj,ki,i], sum(kj), sum(ki)) %*% X[ki,i])

        #logDetV[i] <- with(PCMAbCdEf, det(-2*(A[i,ki,ki]+L[i,ki,ki])))

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
        svdV = svd(matrix(PCMAbCdEf$V[ki,ki,i], sum(ki)), 0, 0)$d
        if(min(svdV)/max(svdV) < threshold_SV) {
          err <- paste0(
            "ERR:02031:PCMBase:MultivariatePCM.R:PCMLmr.default:",i,":",
            " The matrix V for node ", i,
            " is nearly singular: min(svdV)/max(svdV)=", min(svdV)/max(svdV),
            ", det(V)=", det(matrix(PCMAbCdEf$V[i,ki,ki], sum(ki))),
            ". Check the model parameters and the length of the branch",
            " leading to the node. For details on this error, read the User Guide.")
          stop(err)
        }



        # auxilary variables to avoid redundant evaluation
        AplusL <- as.matrix(PCMAbCdEf$A[ki,ki,i] + L[ki,ki,i])
        # ensure symmetry of AplusL, this should guarantee that AplusL_1 is symmetric
        # as well (unless solve-implementation is buggy.)
        AplusL <- 0.5 * (AplusL + t(AplusL))

        AplusL_1 <- solve(AplusL)

        EAplusL_1 <- matrix(PCMAbCdEf$E[kj,ki,i], sum(kj), sum(ki)) %*% AplusL_1
        logDetVNode <- log(det(-2*AplusL))

        # here it is important that we first evaluate r[i] and then m[i,kj]
        # since the expression for r[i] refers to to the value of m[i,ki]
        # before updating it.
        r[i] <- with(PCMAbCdEf, f[i]+r[i]+(sum(ki)/2)*log2pi-.5*logDetVNode -
                       .25*t(b[ki,i]+m[ki,i]) %*% AplusL_1 %*% (b[ki,i]+m[ki,i]))

        m[kj,i] <- with(PCMAbCdEf, d[kj,i] - .5*EAplusL_1 %*% (b[ki,i]+m[ki,i]))

        L[kj,kj,i] <- with(
          PCMAbCdEf,
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
    c(PCMAbCdEf[c("A", "b", "C", "d", "E", "f", "V", "V_1")],
      list(L = L, m = m, r = r))
  }
}



#' Extract error information from a formatted error message.
#' @param x character string representing the error message.
#' @description The function searches x for a pattern matching the format
#' 'ERR:5-digit-code:project-name:source-file:error-specifics:'. Specifically it
#' searches for a regular expression pattern "ERR:[0-9]+:[^:]+:[^:]+:[^:]+:[^:]*:".
#' @return a named list with the parsed error information or NULL, if no match
#' was found. The elements of this list are named as follows:
#' \item{type}{The type of the error message. Usually this is ERROR, but could be
#' WARNING or anything else.}
#' \item{icode}{An integer code of the error.}
#' \item{project}{The name of the project locating the code that raised the error.}
#' \item{file}{The name of the source-file containing the code that raised the error.}
#' \item{fun}{The name of the function raising the error}
#' \item{info}{A character string containing additional error-specific information}
#' \item{msg}{A verbal description of the error.}
#' @export
PCMParseErrorMessage <- function(x) {
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
