source("../PosteriorQuantileCommon.R")

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

# res.validate <- validate.common(
#   modelFromVectorOU, priorOU, genParamOU,
#   n.batch = c(2, 1, 2, 1, 1, 1),
#   params.batch = expression(H["b,ii"],H["b,ij"],theta["b,i"],Sigma["ab,ii"],Sigma["a,12"],Sigma['e,ab,ii']))
#
# save(res.validate, file="BayesValidateOU.RData")



set.seed(1)   # good : 1, 2, 5
param <- genParamOU()
model <- modelFromVectorOU(param)

mvcond(tree, model, r = 1)$vcov(1)
mvcond(tree, model, r = 2)$vcov(1)

X <- mvsim(tree, model, X0=model$X0)
X$values[sample(x=1:length(X$values), 51)] <- NA

gplotX(X, tree) +
  scale_size_continuous(range = c(0.2, 2)) +
  geom_vline(aes(xintercept = model$Theta[1, 1]), col = 2) +
  geom_hline(aes(yintercept = model$Theta[2, 1]), col = 2) +
  geom_vline(aes(xintercept = model$Theta[1, 2]), col = 3) +
  geom_hline(aes(yintercept = model$Theta[2, 2]), col = 3)



