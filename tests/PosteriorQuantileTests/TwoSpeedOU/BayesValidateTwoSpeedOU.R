source("../PosteriorQuantileCommon.R")

modelFromVectorTwoSpeedOU <- function(x){
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
  class(model) <- "TwoSpeedOU"
  model
}

priorTwoSpeedOU <- function(model){
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

genParamTwoSpeedOU <- function() {
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

# res.validate <- validate.common(
#   modelFromVectorTwoSpeedOU, priorTwoSpeedOU, genParamTwoSpeedOU,
#   n.batch = c(1, 1, 1, 2, 1, 2, 1, 1, 1),
#   params.batch = expression(H1["a,ii"],H1["b,ii"],H1["b,ij"],
#                             H2["b,ii"],H2["b,ij"],
#                             theta["b,i"],
#                             Sigma["ab,ii"],Sigma["a,12"],
#                             Sigma['e,ab,ii']))
#
# save(res.validate, file="BayesValidateTwoSpeedOU.RData")



set.seed(5)   # good : 1, 2, 5
param <- genParamTwoSpeedOU()
model <- modelFromVectorTwoSpeedOU(param)

mvcond(tree, model, r = 1)$vcov(1)
mvcond(tree, model, r = 2)$vcov(1)

X <- mvsim(tree, model, X0=model$X0)
X$values[sample(x=1:length(X$values), 51)] <- NA

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


cols <- gg_color_hue(2)
names(cols) <- c("a", "b")

# plotting the tree using the phytools function plotSimmap:
plotSimmap(tree.ab, fsize = 0.01, type="fan", setEnv = TRUE, colors = cols)

gplotX(X, tree) +
  scale_size_continuous(range = c(0.2, 2)) +
  geom_vline(aes(xintercept = model$Theta[1, 1]), col = cols["a"]) +
  geom_hline(aes(yintercept = model$Theta[2, 1]), col = cols["a"]) +
  geom_vline(aes(xintercept = model$Theta[1, 2]), col = cols["b"]) +
  geom_hline(aes(yintercept = model$Theta[2, 2]), col = cols["b"])


