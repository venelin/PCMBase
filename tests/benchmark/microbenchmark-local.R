library(microbenchmark)
library(diversitree)
library(geiger)
library(data.table)
library(phylolm)
library(Rphylopars)
library(PCMBase)
library(PCMBaseCpp)


args <- commandArgs(trailingOnly = TRUE)

if(length(args) > 0) {
  ids <- as.integer(args)
}

compilerInfo <- "icpc-omp-for-simd"
saveResults <- TRUE

if(R.version[['os']]=='linux-gnu') {
  # this only works on linux
  cpuInfo <- system("cat /proc/cpuinfo | grep 'model name' | uniq", intern = TRUE)
} else {
  # this only works on mac OS x
  cpuInfo <- system("sysctl -a -n machdep.cpu.brand_string", intern = TRUE)
}

print(cpuInfo)

nTests <- 5 # number of executions in each microbenchmark

times <- values <- trees <- NULL

# load the trees data.table
load('Trees.RData')

getModeStr <- function(mode) {
  if(mode == 11) "serial" else if(mode == 21) "parallel" else stop("mode not supported");
}

timeRphylopars <- function(X, tree) {
  cat("timeRphylopars... ")
  N <- length(tree$tip.label)
  Z <- X$values + X$errors
  data <- cbind(data.table(species = tree$tip.label), as.data.table(t(Z)[1:N,]))
  time <- system.time(output <- capture.output(
    p_OU <- phylopars(trait_data = data, tree = tree, model="OU", optim_verbose = TRUE)
    )
    )[3]*1000
  #print(output)
  timeLik <- mean(sapply(output, function(o) {
    as.double(strsplit(o, " ")[[1]][3])
  }))
  cat(timeLik, "ms\n")


  if(is.na(timeLik)) {

  }

  timeExtra <-
  data.table(impl = paste0("Rphylopars: C++, 3-point"), time = timeLik, timeExtra = if(is.na(timeLik)) time else time - length(output)*timeLik)
}

timePCMBaseR <- function(X, tree, model) {
  cat("timePCMBaseR... ")
  Z <- X$values + X$errors

  timeExtra <- system.time(
    {
      pruneI <- PCMBase::PCMPruningOrder(tree)
      metaI <- PCMValidate(tree, model)
    })[3]*1000

  timeLik <-
    system.time(
      for(i in 1:nTests)
        PCMBase::PCMLik(Z[, 1:metaI$N,drop=FALSE], tree, model, metaI= metaI, pruneI=pruneI)
)[3]*1000/nTests

  cat(timeLik, "ms\n")

  data.table(impl = paste0("PCMBase: R, Poly"), time = timeLik, timeExtra = timeExtra)
}

timePCMBaseCpp <- function(X, tree, model, mode = 0) {
  cat("timePCMBaseCpp... ")

  options("splittree.postorder.mode" = mode)
  Z <- X$values + X$errors
  timeExtra <- system.time(
    {
      metaI <- PCMValidate(tree, model)
      pruneI <- PCMBaseCpp::PCMCppPruningObject(Z[,1:metaI$N,drop=FALSE], tree, model, metaI)
    })[3]*1000

  timeLik <-
    system.time(
      for(i in 1:nTests)
        PCMBase::PCMLik(Z[,1:metaI$N,drop=FALSE], tree, model, metaI= metaI, pruneI=pruneI)
    )[3]*1000/nTests

  cat(timeLik, "ms\n")

  data.table(impl = paste0("PCMBase: C++, Poly, ", getModeStr(mode)), time = timeLik, timeExtra = timeExtra)
}

timePOUMMR <- function(X, tree, model) {
  cat("timePOUMMR... ")

  N <- length(tree$tip.label)
  Z <- (X$values + X$errors)[,1:N]
  if(is.matrix(Z)) {
    list(impl = paste0("POUMM: R, Poly"), time = as.double(NA), timeExtra = as.double(NA))
    cat(NA, "ms\n")
  } else {
    timeExtra <- system.time(
      {
        pruneI <- POUMM:::pruneTree(tree, Z)
      }
    )[3]*1000

    timeLik <- system.time(
      for(i in 1:(10*nTests)) POUMM::likPOUMMGivenTreeVTips(Z, tree, alpha = model$H[1,1,1], theta=model$Theta[1,1], sigma = model$Sigma[1,1,1], sigmae = model$Sigmae[1,1,1], g0=model$X0[1],  pruneInfo = pruneI)
    )[3]*1000/(10*nTests)
    cat(timeLik, "ms\n")

    data.table(impl = paste0("POUMM: R, Poly"), time = timeLik, timeExtra = timeExtra)
  }
}

timePOUMMCpp <- function(X, tree, model, mode) {
  options("splittree.postorder.mode" = mode)
  N <- length(tree$tip.label)
  Z <- (X$values + X$errors)[,1:N]
  if(is.matrix(Z)) {
    data.table(impl = paste0("POUMM: C++, Poly, ", getModeStr(mode)), time = as.double(NA), timeExtra = as.double(NA))
  } else {
    timeExtra <- system.time(
      {
        pruneI <- POUMM:::pruneTree(tree, Z)
      }
    )[3]*1000

    timeLik <- system.time(
      for(i in 1:(10*nTests)) POUMM::likPOUMMGivenTreeVTipsC(integrator = pruneI$integrator, alpha = model$H[1,1,1], theta=model$Theta[1,1], sigma = model$Sigma[1,1,1], sigmae = model$Sigmae[1,1,1], g0=model$X0[1])
    )[3]*1000/(10*nTests)

    data.table(impl = paste0("POUMM: C++, Poly, ", getModeStr(mode)), time = timeLik, timeExtra = timeExtra)
  }
}

# needed for warmup -loading of C++ libraries
cat("Warming up on id = 1 ")
id <- 1

tree <- trees$tree[[id]]

X <- trees$X[[id]]
N <- trees$N[[id]]
k <- trees$num_traits[[id]]
model <- trees$model[[id]]

print(model)

collessNorm <- trees$collessNorm[[id]]

pSymb <- trees$pSymb[[id]]

x <- timePOUMMR(X, tree, model)
x <- timePCMBaseR(X, tree, model)
x <- timePOUMMCpp(X, tree, model, 11)
x <- timePOUMMCpp(X, tree, model, 21)
x <- timeRphylopars(X=X, tree = tree)
x <- timePCMBaseCpp(X, tree, model, 11)
x <- timePCMBaseCpp(X, tree, model, 21)
rm(x)

gc()

for(id in ids) {
  tree <- trees$tree[[id]]

  X <- trees$X[[id]]
  N <- trees$N[[id]]
  k <- trees$num_traits[[id]]
  model <- trees$model[[id]]

  collessNorm <- trees$collessNorm[[id]]

  pSymb <- trees$pSymb[[id]]

  cat("Tree-id:", id, "; tree-size:", trees$N[[id]], " Colless:", trees$collessNorm[[id]], ", num_traits=", k, "\n")

  times <- rbind(
    times,
    cbind(data.table(compilerInfo = compilerInfo,
                     cpuInfo=cpuInfo, nCores = 4, pSymb = pSymb, collessNorm = collessNorm, N = N, k = k),
          rbindlist(list(
            timePOUMMR(X, tree, model),
            timePCMBaseR(X, tree, model),
            timePOUMMCpp(X, tree, model, 11),
            timePOUMMCpp(X, tree, model, 21),
            timeRphylopars(X=X, tree = tree),
            timePCMBaseCpp(X, tree, model, 11),
            timePCMBaseCpp(X, tree, model, 21)
          ))))

  print(times)
}

if(saveResults) {
  times1 <- copy(times)

  if(!any(ids == 1)) {
    load("Results-microbenchmark-local-1.RData")
  }

  times <- rbind(times, times1)
  save(times,
       file = paste0("Results-microbenchmark-local-1.RData"))
}
