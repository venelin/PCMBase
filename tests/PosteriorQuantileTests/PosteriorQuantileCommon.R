library(adaptMCMC)
library(BayesValidate)
library(PCMBase)
library(PCMBaseCpp)
library(abind)
library(ape)
library(OptimMCMC)
library(coda)

generate.data.common <- function(params, input) {
  model <- input$modelFromVector(params)
  tree <- input$tree
  X <- mvsim(tree, model, X0 = model$X0)
  cat("Generate.data.common called...\n")
  Z <- X$values[, 1:length(tree$tip.label)] + X$errors[, 1:length(tree$tip.label)]

  # some randomly produced list of 51 integers
  NAs <- c(756, 892, 1019, 604, 8, 258, 765, 115, 815, 880, 1013, 313, 788, 794,
           131, 342, 605, 568, 597, 998, 947, 178, 962, 1003, 31,
           702, 52, 344, 296, 925, 211, 963, 678, 142)
  Z[NAs] <- NA
  Z
}

##generate posterior sample
analyze.data.common <- function(X, params.true, input) {
  fit <- runOptimAndMCMC(
    lik = function(par, env.input) {

      pid <- Sys.getpid()
      if(!isTRUE(all.equal(env.input$.RENEW.CACHE, FALSE))) {
        cat("Creating Cpp object on process with pid ", pid, "\n")
        model <- env.input$modelFromVector(par)
        env.input$metaI <- validateModel(env.input$tree, model)
        env.input$pruneInfoObj <- newCppObject(X, env.input$tree, model)
        options("PCMBase.Value.NA"=-1e20)
        options("splittree.postorder.mode" = 11)
        env.input$.RENEW.CACHE <- FALSE
        #save(env.input, file=paste0("env.input.",pid,".RData"))
      }

      value <- mvlik(
            X = NULL,
            tree = NULL,
            model = env.input$modelFromVector(par),
            pruneI = env.input$pruneInfoObj,
            metaI = env.input$metaI,
            log = TRUE)
      #cat("lik(", toString(par), "):", value, "\n")
      value
    },
    prior = function(par, env.input) {
      model <- env.input$modelFromVector(par)
      value <- env.input$prior(model)
      #cat("prior(", toString(par), "):", value, "\n")
      value
    },
    input.data = input,
    config = configOptimAndMCMC(
      parLower = input$parLower,
      parUpper = input$parUpper,
      nCallsOptim = 4,
      listParInitMCMC = lapply(1:4, function(i) input$generate.param()),
      parallelMCMC = FALSE),
    verbose = TRUE)

  save(X, params.true, fit, file = paste0("Replication-ts", as.numeric(Sys.time())*1000, ".RData"))

  lst <- getMCMCs(fit, samplePrior = FALSE)
  as.matrix(window(lst, start = 0.2*end(lst), end=end(lst)))
}

# a common tree for all tests
set.seed(321)

# number of extant species
N=400
# number of regimes
R <- 2
# number of traits
k <- 2

options(digits = 4)

if(require(phytools)) {
  tree.a <- pbtree(b = 1, d = 0.2, n = N)

  # rate mtrix of transition from one regime to another
  Q <- matrix(c(-.1, 0.02, .1, -0.02), R, R)
  colnames(Q) <- rownames(Q) <- letters[1:R]

  tree.ab <- phytools::sim.history(tree.a, Q, anc='a')

  # convert the simmap tree to a normal phylo object with singleton nodes at the
  # within-branch regime changes. The regimes are encoded as names of the
  # edge.length vector
  tree.ab.singles <- map.to.singleton(tree.ab)
  tree <- tree.ab.singles

  tree$edge.regime <- names(tree$edge.length)

  # set jumps at the nodes, where the regime changes
  tree$edge.jump <- rep(0L, length(tree$edge.length))

  # final number of tips (can be different from the previously set value above,
  # because of the non-zero extinction rate in the pbtree call)
  N <- length(tree$tip.label)

  # total number of nodes
  M <- length(tree$edge.length) + 1

  # a loop over all internal nodes
  for(i in (N+2):M) {
    # index of the edge ending at i
    e.parent <- which(tree$edge[,2]==i)

    #cat("i: ", i, "; Parent: ", e.parent, "; ")

    # indices of edges starting at i
    e.daughters <- which(tree$edge[,1]==i)

    #cat(" Children: ", toString(e.daughters), "\n")

    for(j in e.daughters) {
      if(tree$edge.regime[j] != tree$edge.regime[e.parent]) {
        tree$edge.jump[j] <- 1L
      }
    }
  }

  save(tree, file="Tree.RData")
} else {
  load("Tree.RData")
}

# NUMBER OF REPLICATIONS
n.rep <- 48

validate.common <- function(modelFromVector, prior, genParam, ...) {
  par.random <- sapply(1:100, function(i) genParam())

  env.input <- new.env()
  env.input$parLower <- apply(par.random, 1, min)
  env.input$parUpper <- apply(par.random, 1, max)
  env.input$tree <- tree
  env.input$modelFromVector <- modelFromVector
  env.input$prior <- prior
  env.input$generate.param <- genParam


  # set up a parallel cluster on the local computer for parallel MCMC:
  cluster <- parallel::makeCluster(parallel::detectCores(logical = TRUE))
  doParallel::registerDoParallel(cluster)

  res.validate <- validate(
    n.rep = n.rep,
    generate.param = genParam,
    generate.data = generate.data.common,
    generate.data.inputs = env.input,
    analyze.data = analyze.data.common,
    analyze.data.inputs = env.input,
    parallel.rep = TRUE,
    return.all = TRUE,
    ...)

  # Don't forget to destroy the parallel cluster to avoid leaving zombie worker-processes.
  parallel::stopCluster(cluster)
  list(env.input = env.input, res.validate = res.validate)
}


#if(TRUE) {
  # MANUAL EXECUTION ONLY

  library(data.table)
  library(ggplot2)

  gplotX <- function(X, tree) {
    N <- length(tree$tip.label)
    times <- POUMM::nodeTimes(tree, tipsOnly = TRUE)
    Xt <- t(X$values[, 1:N] + X$errors[, 1:N])
    data <- as.data.table(Xt)
    data[, id:=1:(.N)]
    data[, time:=times]
    data[, regime:=sapply(id, function(i) tree$edge.regime[which(tree$edge[, 2]==i)])]

    pl <- ggplot(data) + geom_point(aes(x=V1, y=V2, col=regime, size = time, alpha = time))
    pl
  }

  # plotting the tree using the phytools function plotSimmap:
  plotSimmap(tree.ab, fsize = 0.01, type="fan", setEnv = TRUE)


#}

