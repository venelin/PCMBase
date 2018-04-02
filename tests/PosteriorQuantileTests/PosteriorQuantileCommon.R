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
  X <- PCMSim(tree, model, X0 = model$X0)
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
        env.input$metaI <- PCMValidate(env.input$tree, model)
        env.input$pruneInfoObj <- PCMCppPruningObject(X, env.input$tree, model)
        options("PCMBase.Value.NA"=-1e20)
        options("splittree.postorder.mode" = 11)
        env.input$.RENEW.CACHE <- FALSE
        #save(env.input, file=paste0("env.input.",pid,".RData"))
      }

      value <- PCMLik(
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
      nChainsMCMC = 5,
      samplePriorMCMC = c(TRUE, TRUE, FALSE, FALSE, FALSE),
      listParInitMCMC = lapply(1:5, function(i) input$generate.param()),
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

## Settings for the different models

# BM
modelFromVectorBM <- function(x){
  Sigma.a <- diag(rep(x[1], 2))
  Sigma.a[1,2] <- Sigma.a[2,1] <- x[2]
  Sigma.b <- diag(diag(Sigma.a))

  Sigmae.a <- Sigmae.b <- diag(rep(x[3], 2))

  Sigma <- abind(Sigma.a, Sigma.b, along = 3, new.names = list(x = NULL, y = NULL, regime = c("a", "b")))
  Sigmae <- abind(Sigmae.a, Sigmae.b, along = 3, new.names = list(x = NULL, y = NULL, regime = c("a", "b")))

  model <- list(X0 = c(0, 0), Sigma = Sigma, Sigmae = Sigmae)
  class(model) <- "BM"
  model
}

priorBM <- function(model){
  if(det(model$Sigma[,,1]) < 0 |
     det(model$Sigmae[,,1]) < 0 |
     det(model$Sigma[,,2]) < 0 |
     det(model$Sigmae[,,2]) < 0) {
    -Inf
  } else {
    sum(c(dexp(model$Sigma[1,1,1], 1, log = TRUE),
          dunif(model$Sigma[1,2,1],
                min = -model$Sigma[1,1,1] *.9,
                max = model$Sigma[1,1,1] * .9, log = TRUE),
          dexp(model$Sigmae[1,1,1], 10, log = TRUE)))
  }
}

genParamBM <- function() {
  Sigma.a <- diag(rep(rexp(1, rate = 1), 2))
  Sigma.a[1,2] <- Sigma.a[2,1] <- runif(1, min = -Sigma.a[1,1] * .9, max = Sigma.a[1,1] * .9)
  Sigma.b <- diag(diag(Sigma.a))

  Sigmae.a <- Sigmae.b <- diag(rep(rexp(1, rate = 10), 2))

  c("Sigma[a11]" = Sigma.a[1,1], "Sigma[a12]" = Sigma.a[1,2], "Sigmae[a11]" = Sigmae.a[1,1])
}

# OU
modelFromVectorOU <- function(x){
  H2.a <- matrix(0, 2, 2)
  H2.b <- diag(x[1:2])
  H2.b[1,2] <- H2.b[2,1] <- x[3]

  Theta.a <- c(0, 0)
  Theta.b <- x[4:5]

  Sigma.a <- diag(rep(x[6], 2))
  Sigma.a[1,2] <- Sigma.a[2,1] <- x[7]

  Sigma.b <- diag(diag(Sigma.a))

  Sigmae.a <- Sigmae.b <- diag(rep(x[8], 2))

  H2 <- abind(H2.a, H2.b, along = 3, new.names = list(x = NULL, y = NULL, regime = c("a", "b")))
  Theta <- abind(Theta.a, Theta.b, along = 2, new.names = list(xy = NULL, regime = c("a", "b")))
  Sigma <- abind(Sigma.a, Sigma.b, along = 3, new.names = list(x = NULL, y = NULL, regime = c("a", "b")))
  Sigmae <- abind(Sigmae.a, Sigmae.b, along = 3, new.names = list(x = NULL, y = NULL, regime = c("a", "b")))

  model <- list(X0 = c(0, 0), H = H2, Theta = Theta, Sigma = Sigma, Sigmae = Sigmae)
  class(model) <- "OU"
  model
}

priorOU <- function(model){
  if(det(model$Sigma[,,1]) < 0 |
     det(model$Sigmae[,,1]) < 0 |
     det(model$H[,,2]) < 0 |
     det(model$Sigma[,,2]) < 0 |
     det(model$Sigmae[,,2]) < 0) {
    -Inf
  } else {
    sum(c(dexp(diag(model$H[,,2]), rate = 1, log = TRUE),
          dunif(model$H[1,2,2],
                min = -sqrt(model$H[1,1,2]*model$H[2,2,2])*.9,
                max = sqrt(model$H[1,1,2]*model$H[2,2,2])*.9, log = TRUE),
          dnorm(model$Theta[1,2], 1, .5, log = TRUE),
          dnorm(model$Theta[2,2], 2, .5, log = TRUE),
          dexp(model$Sigma[1,1,1], 1, log = TRUE),
          dunif(model$Sigma[1,2,1],
                min = -model$Sigma[1,1,1] *.9,
                max = model$Sigma[1,1,1] * .9, log = TRUE),
          dexp(model$Sigmae[1,1,1], 10, log = TRUE)))
  }
}

genParamOU <- function() {
  H2.a <- matrix(0, 2, 2)
  H2.b <- diag(rexp(2, rate = 1))
  H2.b[1,2] <- H2.b[2,1] <- runif(1, min = -sqrt(H2.b[1,1]*H2.b[2,2])*.9, max = sqrt(H2.b[1,1]*H2.b[2,2])*.9)

  Theta.b <- c(0, 0) # not part of arguments
  Theta.b <- c(rnorm(1, 1, .5), rnorm(1, 2, .5))

  Sigma.a <- diag(rep(rexp(1, rate = 1), 2))
  Sigma.a[1,2] <- Sigma.a[2,1] <- runif(1, min = -Sigma.a[1,1] * .9, max = Sigma.a[1,1] * .9)
  Sigma.b <- diag(diag(Sigma.a))

  Sigmae.a <- Sigmae.b <- diag(rep(rexp(1, rate = 10), 2))

  c("H2[b11]" = H2.b[1,1], "H2[b22]" = H2.b[2,2], "H2[b12]" = H2.b[1,2],
    "Theta[b1]" = Theta.b[1], "Theta[b2]" = Theta.b[2],
    "Sigma[a11]" = Sigma.a[1,1], "Sigma[a12]" = Sigma.a[1,2],
    "Sigmae[a11]" = Sigmae.a[1,1])
}

# JOU
modelFromVectorJOU <- function(x){
  H2.a <- matrix(0, 2, 2)
  H2.b <- diag(x[1:2])
  H2.b[1,2] <- H2.b[2,1] <- x[3]

  Theta.a <- c(0, 0)
  Theta.b <- x[4:5]

  Sigma.a <- diag(rep(x[6], 2))
  Sigma.a[1,2] <- Sigma.a[2,1] <- x[7]
  Sigma.b <- diag(diag(Sigma.a))

  Sigmae.a <- Sigmae.b <- diag(rep(x[8], 2))

  # mj.a is interpreted as a jump back to regime a, so the
  mj.a <- -Theta.b
  mj.b <- Theta.b

  Sigmaj.a <- diag(1, nrow = 2, ncol = 2)
  Sigmaj.b <- diag(rep(x[9], 2))

  H2 <- abind(H2.a, H2.b, along = 3, new.names = list(x = NULL, y = NULL, regime = c("a", "b")))
  Theta <- abind(Theta.a, Theta.b, along = 2, new.names = list(xy = NULL, regime = c("a", "b")))
  Sigma <- abind(Sigma.a, Sigma.b, along = 3, new.names = list(x = NULL, y = NULL, regime = c("a", "b")))
  Sigmae <- abind(Sigmae.a, Sigmae.b, along = 3, new.names = list(x = NULL, y = NULL, regime = c("a", "b")))
  mj <- abind(mj.a, mj.b, along = 2, new.names = list(xy = NULL, regime = c("a", "b")))
  Sigmaj <- abind(Sigmaj.a, Sigmaj.b, along = 3, new.names = list(x = NULL, y = NULL, regime = c("a", "b")))

  model <- list(X0 = c(0, 0), H = H2, Theta = Theta, Sigma = Sigma, Sigmae = Sigmae, mj = mj, Sigmaj = Sigmaj)
  class(model) <- "JOU"
  model
}

genParamJOU <- function() {
  H2.a <- matrix(0, 2, 2)
  H2.b <- diag(rexp(2, rate = 1))
  H2.b[1,2] <- H2.b[2,1] <- runif(1, min = -sqrt(H2.b[1,1]*H2.b[2,2])*.9, max = sqrt(H2.b[1,1]*H2.b[2,2])*.9)

  Theta.b <- c(0, 0) # not part of arguments
  Theta.b <- c(rnorm(1, 1, .5), rnorm(1, 2, .5))

  Sigma.a <- diag(rep(rexp(1, rate = 1), 2))
  Sigma.a[1,2] <- Sigma.a[2,1] <- runif(1, min = -Sigma.a[1,1] * .9, max = Sigma.a[1,1] * .9)
  Sigma.b <- diag(diag(Sigma.a))

  Sigmae.a <- Sigmae.b <- diag(rep(rexp(1, rate = 10), 2))

  Sigmaj.a <- diag(1, nrow = 2, ncol = 2)
  Sigmaj.b <- diag(rep(rexp(1, rate = 10), 2))


  c("H2[b11]" = H2.b[1,1], "H2[b22]" = H2.b[2,2], "H2[b12]" = H2.b[1,2],
    "Theta[b1]" = Theta.b[1], "Theta[b2]" = Theta.b[2],
    "Sigma[a11]" = Sigma.a[1,1], "Sigma[a12]" = Sigma.a[1,2],
    "Sigmae[a11]" = Sigmae.a[1,1],
    "Sigmaj[b11]" = Sigmaj.b[1,1])
}

priorJOU <- function(model){
  if(det(model$Sigma[,,1]) < 0 |
     det(model$Sigmae[,,1]) < 0 |
     det(model$Sigmaj[,,1]) < 0 |
     det(model$H[,,2]) < 0 |
     det(model$Sigma[,,2]) < 0 |
     det(model$Sigmae[,,2]) < 0 |
     det(model$Sigmaj[,,2]) < 0) {
    -Inf
  } else {
    sum(c(dexp(diag(model$H[,,2]), rate = 1, log = TRUE),
          dunif(model$H[1,2,2],
                min = -sqrt(model$H[1,1,2] * model$H[2,2,2]) * .9,
                max = sqrt(model$H[1,1,2] * model$H[2,2,2]) * .9, log = TRUE),
          dnorm(model$Theta[1,2], 1, .5, log = TRUE),
          dnorm(model$Theta[2,2], 2, .5, log = TRUE),
          dexp(model$Sigma[1,1,1], 1, log = TRUE),
          dunif(model$Sigma[1,2,1],
                min = -model$Sigma[1,1,1] * .9,
                max = model$Sigma[1,1,1] * .9, log = TRUE),
          dexp(model$Sigmae[1,1,1], 10, log = TRUE),
          dexp(model$Sigmaj[1,1,2], 10, log = TRUE)))
  }
}


# DOU
modelFromVectorDOU <- function(x){
  H1.a <- diag(x[1], nrow = 2, ncol = 2)
  H1.b <- diag(x[2], nrow = 2, ncol = 2)
  H1.b[1,2] <- H1.b[2,1] <- x[3]

  H2.a <- matrix(0, 2, 2)
  H2.b <- diag(x[4:5])
  H2.b[1,2] <- H2.b[2,1] <- x[6]

  Theta.a <- c(0, 0)
  Theta.b <- x[7:8]

  Sigma.a <- diag(rep(x[9], 2))
  Sigma.a[1,2] <- Sigma.a[2,1] <- x[10]

  Sigma.b <- diag(diag(Sigma.a))

  Sigmae.a <- Sigmae.b <- diag(rep(x[11], 2))

  H1 <- abind(H1.a, H1.b, along = 3, new.names = list(x = NULL, y = NULL, regime = c("a", "b")))
  H2 <- abind(H2.a, H2.b, along = 3, new.names = list(x = NULL, y = NULL, regime = c("a", "b")))
  Theta <- abind(Theta.a, Theta.b, along = 2, new.names = list(xy = NULL, regime = c("a", "b")))
  Sigma <- abind(Sigma.a, Sigma.b, along = 3, new.names = list(x = NULL, y = NULL, regime = c("a", "b")))
  Sigmae <- abind(Sigmae.a, Sigmae.b, along = 3, new.names = list(x = NULL, y = NULL, regime = c("a", "b")))

  model <- list(X0 = c(0, 0), H1 = H1, H2 = H2, Theta = Theta, Sigma = Sigma, Sigmae = Sigmae)
  class(model) <- "DOU"
  model
}

priorDOU <- function(model){
  if(det(model$Sigma[,,1]) < 0 |
     det(model$Sigmae[,,1]) < 0 |
     det(model$H2[,,2]) < 0 |
     det(model$Sigma[,,2]) < 0 |
     det(model$Sigmae[,,2]) < 0) {
    -Inf
  } else {
    sum(c(dnorm(model$H1[1,1,1], 0, .5, log = TRUE),
          dnorm(model$H1[1,1,2], 0, .5, log = TRUE),
          dnorm(model$H1[1,2,2], 0, .5, log = TRUE),
          dexp(diag(model$H2[,,2]), rate = 1, log = TRUE),
          dunif(model$H2[1,2,2],
                min = -sqrt(model$H2[1,1,2] * model$H2[2,2,2]) * .9,
                max = sqrt(model$H2[1,1,2] * model$H2[2,2,2]) * .9, log = TRUE),
          dnorm(model$Theta[1,2], 1, .5, log = TRUE),
          dnorm(model$Theta[2,2], 2, .5, log = TRUE),
          dexp(model$Sigma[1,1,1], 1, log = TRUE),
          dunif(model$Sigma[1,2,1],
                min = -model$Sigma[1,1,1] * .9,
                max = model$Sigma[1,1,1] * .9, log = TRUE),
          dexp(model$Sigmae[1,1,1], 10, log = TRUE)))
  }
}

genParamDOU <- function() {
  H1.a <- diag(rnorm(1, 0, .5), nrow = 2, ncol = 2)
  H1.b <- diag(rnorm(1, 0, .5), nrow = 2, ncol = 2)
  H1.b[1,2] <- H1.b[2,1] <- rnorm(1, 0, .5)

  H2.a <- matrix(0, 2, 2)
  H2.b <- diag(rexp(2, rate = 1))
  H2.b[1,2] <- H2.b[2,1] <- runif(1, min = -sqrt(H2.b[1,1]*H2.b[2,2])*.9, max = sqrt(H2.b[1,1]*H2.b[2,2])*.9)

  Theta.b <- c(0, 0) # not part of arguments
  Theta.b <- c(rnorm(1, 1, .5), rnorm(1, 2, .5))

  Sigma.a <- diag(rep(rexp(1, rate = 1), 2))
  Sigma.a[1,2] <- Sigma.a[2,1] <- runif(1, min = -Sigma.a[1,1] * .9, max = Sigma.a[1,1] * .9)
  Sigma.b <- diag(diag(Sigma.a))

  Sigmae.a <- Sigmae.b <- diag(rep(rexp(1, rate = 10), 2))

  c("H1[a11]" = H1.a[1,1], "H1[b11]" = H1.b[1,1], "H1[b12]" = H1.b[1,2],
    "H2[b11]" = H2.b[1,1], "H2[b22]" = H2.b[2,2], "H2[b12]" = H2.b[1,2],
    "Theta[b1]" = Theta.b[1], "Theta[b2]" = Theta.b[2],
    "Sigma[a11]" = Sigma.a[1,1], "Sigma[a12]" = Sigma.a[1,2],
    "Sigmae[a11]" = Sigmae.a[1,1])
}


