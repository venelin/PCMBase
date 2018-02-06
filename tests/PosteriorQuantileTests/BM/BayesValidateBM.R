source("../PosteriorQuantileCommon.R")

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

res.validate <- validate.common(
  modelFromVectorBM, priorBM, genParamBM,
  n.batch = c(1, 1, 1),
  params.batch = expression(Sigma["ab,ii"],Sigma["a,12"],Sigma['e,ab,ii']))

save(res.validate, file="BayesValidateBM.RData")



# set.seed(5)   # good : 1, 2, 5
# param <- genParamOU()
# model <- modelFromVectorOU(param)
#
# mvcond(tree, model, r = 1)$vcov(1)
# mvcond(tree, model, r = 2)$vcov(1)
#
# X <- mvsim(tree, model, X0=model$X0)
# X$values[sample(x=1:length(X$values), 51)] <- NA
#
# plotX(X, tree)
