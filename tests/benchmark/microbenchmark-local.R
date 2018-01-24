library(microbenchmark)
library(diversitree)
library(geiger)
library(data.table)
library(phylolm)
library(SPLiTTree)


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

set.seed(1)
nTests <- 100 # number of executions in each microbenchmark

g0 <- 16
alpha <- 2
theta <- 4
sigma <- .2
sigmae <- .8

a <- rexp(nTests)
th <- runif(nTests, 0, 4)
s <- rexp(nTests)
se <- rexp(nTests)
se[1:(nTests / 2)] <- sigmae

times <- values <- trees <- NULL

# load the trees data.table
load('Trees.RData')

for(id in ids) {
  treeNo <- id
  
  cat("Tree-id:", id, "; tree-size:", trees$N[[id]], " Colless:", trees$collessNorm[[id]], "\n")
  
  tree <- trees$tree[[id]]
  
  z <- trees$z[[id]]
  N <- trees$N[[id]]
  
  collessNorm <- trees$collessNorm[[id]]
  
  th[1:(nTests / 2)] <- mean(z[1:N])
  
  type <- trees$treeType[[id]]
  pSymb <- trees$pSymb[[id]]
  
  time.create.geiger <- time.create.divtree.R <- time.create.divtree.C <- 
    time.create.pruneInfo <- time.create.POUMM_lnDetV_Q <- NA
  
  cat("Creating cache structure for geiger...\n")
  time.create.geiger <- system.time(
    lik.geiger.C <- try({
      geiger:::bm.lik(tree, z[1:N], SE = NA, model="OU")
      },
                        silent = TRUE))
  
  
  cat("Creating cache structure for diversitree R...\n")
  time.create.divtree.R <- system.time(
    lik.divtree.R <- 
      try({
        if(N==1e5 & collessNorm>0.9) {
          stop("stack overflow error")
        }
        make.ou(tree, z[1:N], states.sd =  sigmae, with.optimum = TRUE, 
                control = list(method = 'pruning', backend = 'R'))
        },
          silent = TRUE))
  
  cat("Creating cache structure for diversitree C...\n")
  time.create.divtree.C <- system.time(
    lik.divtree.C <- 
      try({
        if(N==1e5 & collessNorm>0.9) {
          stop("stack overflow error")
        }
        make.ou(tree, z[1:N], states.sd =  sigmae, with.optimum = TRUE, 
                  control = list(method = 'pruning', backend = 'C'))
        },
          silent = TRUE))
  
  cat("Creating pruning object SPLiTTree:::POUMM_abc...\n")
  time.create.POUMM_abc <- system.time(
    POUMM_abc <- SPLiTTree:::ParallelPruningAbcPOUMM$new(tree, z[1:N], rep(0, N)))
  
  cat("Creating pruning object SPLiTTree:::POUMM_lnDetV_Q_1d...\n")
  time.create.POUMM_lnDetV_Q <- system.time(
    POUMM_lnDetV_Q <- SPLiTTree:::ParallelPruningThreePointPOUMM$new(tree, z[1:N], rep(0, N)))
  
  DoPruning_POUMM_abc <- POUMM_abc$TraverseTree
  lik_POUMM_abc <- function(
    objPOUMM_abc, g0, alpha, theta, sigma, sigmae, mode) {
    
    abc<-DoPruning_POUMM_abc(c(alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae), mode)
    
    g0_theta <- g0 - theta
    (abc[1] * g0_theta + abc[2]) * g0_theta + abc[3]
  }
  
  LOG_2PI <- log(2*pi)
  DoPruning_POUMM_lnDetV_Q <- POUMM_lnDetV_Q$TraverseTree
  lik_POUMM_lnDetV_Q <- function(
    objPOUMM_lnDetV_Q, g0, alpha, theta, sigma, sigmae, mode) {
    
    res <- DoPruning_POUMM_lnDetV_Q(
      c(g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae),
      mode)
    
    -1/2*(N * LOG_2PI + 2*res[3] + res[1] + res[2])
  }
  
  anc <- rep(0.0, nTests)
  res_geiger_C <- res_divtree_R <- res_divtree_C <- 
    res_POUMM_abc_0 <- res_POUMM_abc_1 <- res_POUMM_abc_2 <- res_POUMM_abc_3 <- res_POUMM_abc_4 <- res_POUMM_abc_5 <- res_POUMM_abc_6 <- res_POUMM_abc_7 <- res_POUMM_abc_8 <- res_POUMM_abc_9 <- res_POUMM_abc_10 <-
    res_POUMM_lnDetV_Q_0 <- res_POUMM_lnDetV_Q_1 <- res_POUMM_lnDetV_Q_2 <- res_POUMM_lnDetV_Q_3 <- res_POUMM_lnDetV_Q_4 <- res_POUMM_lnDetV_Q_5 <- res_POUMM_lnDetV_Q_6 <- res_POUMM_lnDetV_Q_7 <- res_POUMM_lnDetV_Q_8 <- res_POUMM_lnDetV_Q_9 <- res_POUMM_lnDetV_Q_10 <-
    res_geiger_C <- rep(0.0, nTests)
  
  counter_geiger_C <- 
    counter_divtree_R <- counter_divtree_C <- 
    counter_POUMM_abc_0 <- counter_POUMM_abc_1 <- counter_POUMM_abc_2 <- counter_POUMM_abc_3 <- counter_POUMM_abc_4 <- counter_POUMM_abc_5 <- counter_POUMM_abc_6 <- counter_POUMM_abc_7 <- counter_POUMM_abc_8 <- counter_POUMM_abc_9 <- counter_POUMM_abc_10 <-
    counter_POUMM_lnDetV_Q_0 <- counter_POUMM_lnDetV_Q_1 <- counter_POUMM_lnDetV_Q_2 <- counter_POUMM_lnDetV_Q_3 <- counter_POUMM_lnDetV_Q_4 <- counter_POUMM_lnDetV_Q_5 <- counter_POUMM_lnDetV_Q_6 <- counter_POUMM_lnDetV_Q_7 <- counter_POUMM_lnDetV_Q_8 <- counter_POUMM_lnDetV_Q_9 <- counter_POUMM_lnDetV_Q_10 <- 1
  
  counter_gc <- 1;
  
  
  if(!is.null(lik.geiger.C) & class(lik.geiger.C) != 'try-error') {
    mb_geiger <- microbenchmark(
      "gc" = {if(counter_gc %% 10 == 0) {gc()}; counter_gc <-  counter_gc + 1; },
      "geiger: C++" = {i <- counter_geiger_C; counter_geiger_C <- counter_geiger_C + 1; res_geiger_C[i] <- lik.geiger.C(c(alpha = a[i], sigsq = s[i]^2, SE = se[i])); },
      times = nTests)
  } else {
    cat("lik.geiger not done\n")
    print(lik.geiger.C)
    mb_geiger <- NULL
  }
  
  if(!is.null(lik.divtree.R) & class(lik.divtree.R) != 'try-error') {
    mb_diversitree_R <- microbenchmark(
      "gc" = {if(counter_gc %% 10 == 0) {gc()}; counter_gc <-  counter_gc + 1; },
      "diversitree: R" = {i <- counter_divtree_R;counter_divtree_R <- counter_divtree_R + 1;res_divtree_R[i] <- lik.divtree.R(c(s2=s[i]^2, alpha=a[i], theta=th[i]));},
      times = nTests)
  } else {
    cat("lik.divtree.R not done\n")
    print(lik.divtree.R)
    mb_diversitree_R <- NULL
  }
  
  if(!is.null(lik.divtree.C) & class(lik.divtree.C) != 'try-error') {
    mb_diversitree_C <- microbenchmark("gc" = {if(counter_gc %% 10 == 0) {gc()}; counter_gc <-  counter_gc + 1;},
      "diversitree: C++" = {i <- counter_divtree_C; counter_divtree_C <- counter_divtree_C + 1; res_divtree_C[i] <- lik.divtree.C(c(s2=s[i]^2, alpha=a[i], theta=th[i])); },
      times = nTests)
  } else {
    cat("lik.divtree.C not done\n")
    print(lik.divtree.C)
    mb_diversitree_C <- NULL
  }

  counter_gc <- 1;
  
  mbParallel <- microbenchmark(
    "gc" = { if(counter_gc %% 4 == 0) {gc()}; counter_gc <- counter_gc + 1;},
    "ppa: C++, abc, auto" = {i <- counter_POUMM_abc_0; counter_POUMM_abc_0 <- counter_POUMM_abc_0 + 1; res_POUMM_abc_0[i] <- lik_POUMM_abc(POUMM_abc, g0 = g0, alpha = a[i], theta = th[i],sigma = s[i], sigmae = se[i], mode = 0); },
    "ppa: C++, abc, serial postorder" = {i <- counter_POUMM_abc_1; counter_POUMM_abc_1 <- counter_POUMM_abc_1 + 1; res_POUMM_abc_1[i] <- lik_POUMM_abc(POUMM_abc, g0 = g0, alpha = a[i], theta = th[i],sigma = s[i], sigmae = se[i], mode = 10); },
    "ppa: C++, abc, serial prune-ranges" = {i <- counter_POUMM_abc_2; counter_POUMM_abc_2 <- counter_POUMM_abc_2 + 1; res_POUMM_abc_2[i] <- lik_POUMM_abc(POUMM_abc, g0 = g0, alpha = a[i], theta = th[i],sigma = s[i], sigmae = se[i], mode = 11); },
    "ppa: C++, abc, parallel prune-ranges" = {i <- counter_POUMM_abc_3; counter_POUMM_abc_3 <- counter_POUMM_abc_3 + 1; res_POUMM_abc_3[i] <- lik_POUMM_abc(POUMM_abc, g0 = g0, alpha = a[i], theta = th[i],sigma = s[i], sigmae = se[i], mode = 21); },
    "ppa: C++, abc, parallel split" = {i <- counter_POUMM_abc_4; counter_POUMM_abc_4 <- counter_POUMM_abc_4 + 1; res_POUMM_abc_4[i] <- lik_POUMM_abc(POUMM_abc, g0 = g0, alpha = a[i], theta = th[i],sigma = s[i], sigmae = se[i], mode = 23); },
    "ppa: C++, abc, parallel visit-queue" = {i <- counter_POUMM_abc_5; counter_POUMM_abc_5 <- counter_POUMM_abc_5 + 1; res_POUMM_abc_5[i] <- lik_POUMM_abc(POUMM_abc, g0 = g0, alpha = a[i], theta = th[i],sigma = s[i], sigmae = se[i], mode = 24); },
    "ppa: C++, abc, parallel visit-ranges" = {i <- counter_POUMM_abc_6; counter_POUMM_abc_6 <- counter_POUMM_abc_6 + 1; res_POUMM_abc_6[i] <- lik_POUMM_abc(POUMM_abc, g0 = g0, alpha = a[i], theta = th[i],sigma = s[i], sigmae = se[i], mode = 22); },
    "ppa: C++, abc, hybrid prune-ranges" = {i <- counter_POUMM_abc_7; counter_POUMM_abc_7 <- counter_POUMM_abc_7 + 1; res_POUMM_abc_7[i] <- lik_POUMM_abc(POUMM_abc, g0 = g0, alpha = a[i], theta = th[i],sigma = s[i], sigmae = se[i], mode = 31); },
    "ppa: C++, abc, hybrid split" = {i <- counter_POUMM_abc_8; counter_POUMM_abc_8 <- counter_POUMM_abc_8 + 1; res_POUMM_abc_8[i] <- lik_POUMM_abc(POUMM_abc, g0 = g0, alpha = a[i], theta = th[i],sigma = s[i], sigmae = se[i], mode = 33); },
    "ppa: C++, abc, hybrid visit-ranges" = {i <- counter_POUMM_abc_9; counter_POUMM_abc_9 <- counter_POUMM_abc_9 + 1; res_POUMM_abc_9[i] <- lik_POUMM_abc(POUMM_abc, g0 = g0, alpha = a[i], theta = th[i],sigma = s[i], sigmae = se[i], mode = 32); },
    "ppa: C++, abc, serial visit-ranges" = {i <- counter_POUMM_abc_10; counter_POUMM_abc_10 <- counter_POUMM_abc_10 + 1; res_POUMM_abc_10[i] <- lik_POUMM_abc(POUMM_abc, g0 = g0, alpha = a[i], theta = th[i],sigma = s[i], sigmae = se[i], mode = 12); },
    
    "ppa: C++, lnDetV_Q, auto" = {i <- counter_POUMM_lnDetV_Q_0; counter_POUMM_lnDetV_Q_0 <- counter_POUMM_lnDetV_Q_0 + 1; res_POUMM_lnDetV_Q_0[i] <- lik_POUMM_lnDetV_Q(POUMM_lnDetV_Q, g0 = g0, alpha = a[i], theta = th[i],sigma = s[i], sigmae = se[i], mode = 0); },
    "ppa: C++, lnDetV_Q, serial postorder" = {i <- counter_POUMM_lnDetV_Q_1; counter_POUMM_lnDetV_Q_1 <- counter_POUMM_lnDetV_Q_1 + 1; res_POUMM_lnDetV_Q_1[i] <- lik_POUMM_lnDetV_Q(POUMM_lnDetV_Q, g0 = g0, alpha = a[i], theta = th[i],sigma = s[i], sigmae = se[i], mode = 10); },
    "ppa: C++, lnDetV_Q, serial prune-ranges" = {i <- counter_POUMM_lnDetV_Q_2; counter_POUMM_lnDetV_Q_2 <- counter_POUMM_lnDetV_Q_2 + 1; res_POUMM_lnDetV_Q_2[i] <- lik_POUMM_lnDetV_Q(POUMM_lnDetV_Q, g0 = g0, alpha = a[i], theta = th[i],sigma = s[i], sigmae = se[i], mode = 11); },
    "ppa: C++, lnDetV_Q, parallel prune-ranges" = {i <- counter_POUMM_lnDetV_Q_3; counter_POUMM_lnDetV_Q_3 <- counter_POUMM_lnDetV_Q_3 + 1; res_POUMM_lnDetV_Q_3[i] <- lik_POUMM_lnDetV_Q(POUMM_lnDetV_Q, g0 = g0, alpha = a[i], theta = th[i],sigma = s[i], sigmae = se[i], mode = 21); },
    "ppa: C++, lnDetV_Q, parallel split" = {i <- counter_POUMM_lnDetV_Q_4; counter_POUMM_lnDetV_Q_4 <- counter_POUMM_lnDetV_Q_4 + 1; res_POUMM_lnDetV_Q_4[i] <- lik_POUMM_lnDetV_Q(POUMM_lnDetV_Q, g0 = g0, alpha = a[i], theta = th[i],sigma = s[i], sigmae = se[i], mode = 23); },
    "ppa: C++, lnDetV_Q, parallel visit-queue" = {i <- counter_POUMM_lnDetV_Q_5; counter_POUMM_lnDetV_Q_5 <- counter_POUMM_lnDetV_Q_5 + 1; res_POUMM_lnDetV_Q_5[i] <- lik_POUMM_lnDetV_Q(POUMM_lnDetV_Q, g0 = g0, alpha = a[i], theta = th[i],sigma = s[i], sigmae = se[i], mode = 24); },
    "ppa: C++, lnDetV_Q, parallel visit-ranges" = {i <- counter_POUMM_lnDetV_Q_6; counter_POUMM_lnDetV_Q_6 <- counter_POUMM_lnDetV_Q_6 + 1; res_POUMM_lnDetV_Q_6[i] <- lik_POUMM_lnDetV_Q(POUMM_lnDetV_Q, g0 = g0, alpha = a[i], theta = th[i],sigma = s[i], sigmae = se[i], mode = 22); },
    "ppa: C++, lnDetV_Q, hybrid prune-ranges" = {i <- counter_POUMM_lnDetV_Q_7; counter_POUMM_lnDetV_Q_7 <- counter_POUMM_lnDetV_Q_7 + 1; res_POUMM_lnDetV_Q_7[i] <- lik_POUMM_lnDetV_Q(POUMM_lnDetV_Q, g0 = g0, alpha = a[i], theta = th[i],sigma = s[i], sigmae = se[i], mode = 31); },
    "ppa: C++, lnDetV_Q, hybrid split" = {i <- counter_POUMM_lnDetV_Q_8; counter_POUMM_lnDetV_Q_8 <- counter_POUMM_lnDetV_Q_8 + 1; res_POUMM_lnDetV_Q_8[i] <- lik_POUMM_lnDetV_Q(POUMM_lnDetV_Q, g0 = g0, alpha = a[i], theta = th[i],sigma = s[i], sigmae = se[i], mode = 33); },
    "ppa: C++, lnDetV_Q, hybrid visit-ranges" = {i <- counter_POUMM_lnDetV_Q_9; counter_POUMM_lnDetV_Q_9 <- counter_POUMM_lnDetV_Q_9 + 1; res_POUMM_lnDetV_Q_9[i] <- lik_POUMM_lnDetV_Q(POUMM_lnDetV_Q, g0 = g0, alpha = a[i], theta = th[i],sigma = s[i], sigmae = se[i], mode = 32); },
    "ppa: C++, lnDetV_Q, serial visit-ranges" = {i <- counter_POUMM_lnDetV_Q_10; counter_POUMM_lnDetV_Q_10 <- counter_POUMM_lnDetV_Q_10 + 1; res_POUMM_lnDetV_Q_10[i] <- lik_POUMM_lnDetV_Q(POUMM_lnDetV_Q, g0 = g0, alpha = a[i], theta = th[i],sigma = s[i], sigmae = se[i], mode = 12); },
    times = nTests
  )
  
  values <- 
    rbind(values, 
          data.table(compilerInfo = compilerInfo, 
                     cpuInfo=cpuInfo, nCores = NA,
                     type = type, pSymb = pSymb, collessNorm = collessNorm,  
                     treeNo = treeNo, N = N, 
                     testId = 1:nTests, 
                     alpha = a, theta = th, sigma = s, sigmae = se, 
                     "geiger: C++" = res_geiger_C, 
                     "diversitree: R" = res_divtree_R, 
                     "diversitree: C++" = res_divtree_C,  
                     
                     "ppa: C++, abc, auto" = res_POUMM_abc_0,
                     "ppa: C++, abc, serial postorder" = res_POUMM_abc_1,
                     "ppa: C++, abc, serial prune-ranges" = res_POUMM_abc_2,
                     "ppa: C++, abc, parallel prune-ranges" = res_POUMM_abc_3,
                     "ppa: C++, abc, parallel split" = res_POUMM_abc_4,
                     "ppa: C++, abc, parallel visit-queue" = res_POUMM_abc_5,
                     "ppa: C++, abc, parallel visit-ranges" = res_POUMM_abc_6,
                     "ppa: C++, abc, hybrid prune-ranges" = res_POUMM_abc_7,
                     "ppa: C++, abc, hybrid split" = res_POUMM_abc_8,
                     "ppa: C++, abc, hybrid visit-ranges" = res_POUMM_abc_9,
                     "ppa: C++, abc, serial visit-ranges" = res_POUMM_abc_10,
                     
                     "ppa: C++, lnDetV_Q, auto" = res_POUMM_lnDetV_Q_0,
                     "ppa: C++, lnDetV_Q, serial postorder" = res_POUMM_lnDetV_Q_1,
                     "ppa: C++, lnDetV_Q, serial prune-ranges" = res_POUMM_lnDetV_Q_2,
                     "ppa: C++, lnDetV_Q, parallel prune-ranges" = res_POUMM_lnDetV_Q_3,
                     "ppa: C++, lnDetV_Q, parallel split" = res_POUMM_lnDetV_Q_4,
                     "ppa: C++, lnDetV_Q, parallel visit-queue" = res_POUMM_lnDetV_Q_5,
                     "ppa: C++, lnDetV_Q, parallel visit-ranges" = res_POUMM_lnDetV_Q_6,
                     "ppa: C++, lnDetV_Q, hybrid prune-ranges" = res_POUMM_lnDetV_Q_7,
                     "ppa: C++, lnDetV_Q, hybrid split" = res_POUMM_lnDetV_Q_8,
                     "ppa: C++, lnDetV_Q, hybrid visit-ranges" = res_POUMM_lnDetV_Q_9,
                     "ppa: C++, lnDetV_Q, serial visit-ranges" = res_POUMM_lnDetV_Q_10))
  
    times <- rbind(times, 
          cbind(data.table(compilerInfo = compilerInfo, 
                           cpuInfo=cpuInfo, nCores = 4, 
                           type = type, pSymb = pSymb, collessNorm = collessNorm,  
                           N = N, 
                           timeCreateCache = c(NA, 
                                               rep(time.create.POUMM_abc[3], 10),
                                               rep(time.create.POUMM_lnDetV_Q[3], 10))),
                           as.data.table(print(mbParallel, unit="ms"))))
  
  
  if(!is.null(mb_geiger) & class(mb_geiger) != "try-error") {
    times <- 
      rbind(times, 
            cbind(data.table(compilerInfo = compilerInfo, 
                             cpuInfo=cpuInfo, nCores = 1, 
                             type = type, pSymb = pSymb, collessNorm = collessNorm,  
                             N = N, 
                             timeCreateCache = c(NA, time.create.geiger[3])),
                             as.data.table(print(mb_geiger, unit="ms"))))
  }
  if(!is.null(mb_diversitree_R) & class(mb_diversitree_R) != "try-error") {
    times <- 
      rbind(times, 
            cbind(data.table(compilerInfo = compilerInfo, 
                             cpuInfo=cpuInfo, nCores = 1, 
                             type = type, pSymb = pSymb, collessNorm = collessNorm,  
                             N = N, 
                             timeCreateCache = c(NA, time.create.divtree.R[3])), 
                             as.data.table(print(mb_diversitree_R, unit="ms"))))
  }
  if(!is.null(mb_diversitree_C) & class(mb_diversitree_C) != "try-error") {
    times <- 
      rbind(times, 
            cbind(data.table(compilerInfo = compilerInfo, 
                             cpuInfo=cpuInfo, nCores = 1, 
                             type = type, pSymb = pSymb, collessNorm = collessNorm,  
                             N = N, 
                             timeCreateCache = c(NA, time.create.divtree.C[3])),
                             as.data.table(print(mb_diversitree_C, unit="ms"))))
  }
}

if(saveResults) {
  values1 <- copy(values)
  times1 <- copy(times)

  if(!any(ids == 1)) {
    load("Results-microbenchmark-local-6.RData")
  }
  
  values <- rbind(values, values1)
  times <- rbind(times, times1)
  save(values, times, 
       file = paste0("Results-microbenchmark-local-6.RData"))
}
