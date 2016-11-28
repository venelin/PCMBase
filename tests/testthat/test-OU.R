library(phytools)
library(testthat)


set.seed(1)


#numerical zero:
EPS <- 10^-9

# number of regimes
R <- 2
# number of traits
k <- 3


# rate mtrix of transition from one regime to another
Q <- matrix(c(-1, 1, 1, -1), R, R)
colnames(Q) <- rownames(Q) <- letters[1:R]

## Specifying a bivariate OU process for each regime
# First, specify the A, theta, Sigma and sigmae2 parameters for each regime.
# Then we use the abind function to stack the parameters into arrays which's first
# dimension is the regime

# regimes

# in regime 'a' the two traits evolve according to two independent OU processes
a.A <- rbind(c(2, 0, 0),
             c(0, .6, 0),
             c(0, 0, .2))
a.theta <- c(10, 6, 2)
a.Sigma <- rbind(c(1.6, 0, 0),
                 c(0, 2.4, 0),
                 c(0, 0, 2))
a.Sigmae2 <- rbind(c(0, 0, 0),
                   c(0, 0, 0),
                   c(0, 0, 0))

# in regime 'b' there is correlation between the traits
b.A <- rbind(c(2, .1, .2),
             c(.1, .6, .2),
             c(.2, .2, .3))
b.theta <- c(10, 6, 2)
b.Sigma <- rbind(c(1.6, .3, .3),
                 c(.3, 0.3, .4),
                 c(.3, .4, 2))
b.Sigmae2 <- rbind(c(.2, 0, 0),
                   c(0, .3, 0),
                   c(0, 0, .4))

A <- abind::abind(a.A, b.A, along=-1, new.names=list(regime=c('a','b'), x=NULL, y=NULL))
Theta <- abind::abind(a.theta, b.theta, along=-1, new.names=list(regime=c('a', 'b'), xy=NULL))
Sigma <- abind::abind(a.Sigma, b.Sigma, along=-1, new.names=list(regime=c('a','b'), x=NULL, y=NULL))
Sigmae <- abind::abind(a.Sigmae2, b.Sigmae2, along=-1, new.names=list(regime=c('a','b'), x=NULL, y=NULL))


## Simulations of trait data

# regime 'a', trait 1
params.a.1 <- list(A=A['a',1,1,drop=FALSE],
Theta=Theta['a',1,drop=FALSE],
Sigma=Sigma['a',1,1,drop=FALSE],
Sigmae=Sigmae['a',1,1,drop=FALSE])
class(params.a.1) <- 'OU'

# regime 'a', trait 2
params.a.2 <- list(A=A['a',2,2,drop=FALSE],
Theta=Theta['a',2,drop=FALSE],
Sigma=Sigma['a',2,2,drop=FALSE],
Sigmae=Sigmae['a',2,2,drop=FALSE])
class(params.a.2) <- 'OU'

# regime 'a', trait 3
params.a.3 <- list(A=A['a',3,3,drop=FALSE],
                   Theta=Theta['a',3,drop=FALSE],
                   Sigma=Sigma['a',3,3,drop=FALSE],
                   Sigmae=Sigmae['a',3,3,drop=FALSE])
class(params.a.3) <- 'OU'

# regime 'a', traits 1, 2 and 3
params.a.123 <- list(A=A['a',,,drop=FALSE],
Theta=Theta['a',,drop=FALSE],
Sigma=Sigma['a',,,drop=FALSE],
Sigmae=Sigmae['a',,,drop=FALSE])
class(params.a.123) <- 'OU'


# regime 'b', traits 1 and 2
params.b.123 <- list(A=A['b',,,drop=FALSE],
Theta=Theta['b',,drop=FALSE],
Sigma=Sigma['b',,,drop=FALSE],
Sigmae=Sigmae['b',,,drop=FALSE])
class(params.b.123) <- 'OU'

# regimes 'a' and 'b', traits 1 and 2
params.ab.123 <- list(A=A[,,,drop=FALSE],
Theta=Theta[,,drop=FALSE],
Sigma=Sigma[,,,drop=FALSE],
Sigmae=Sigmae[,,,drop=FALSE])
class(params.ab.123) <- 'OU'



context(ctx <- "R=1/k=1/N=2")

# number of tips
N <- 2

# tree with one regime
tree.a <- phytools::pbtree(n=N, scale=1)
# tree with two regimes

params.a.1$tree <- params.a.2$tree <- params.a.3$tree <- tree.a

# generate traits
traits.a.1 <- generateMVTrait(params.a.1, 0, verbose=TRUE)
traits.a.2 <- generateMVTrait(params.a.2, 0, verbose=TRUE)
traits.a.3 <- generateMVTrait(params.a.2, 0, verbose=TRUE)

# test likelihood
test_that(paste(ctx, "Match univariate likelihood from patherit regime a"), {
  expect_true(
    abs(lik(traits.a.1$values+traits.a.1$errors,
            params.a.1) -
          patherit::lik.poumm(traits.a.1$values[,1]+traits.a.1$errors[,1],
                              params.a.1$tree,
                              params.a.1$A[1,1,1],
                              params.a.1$Theta[1,1],
                              sqrt(params.a.1$Sigma[1,1,1]),
                              sqrt(params.a.1$Sigmae[1,1,1]),
                              distgr='maxlik', usempfr=0)) < EPS)
  expect_true(abs(
    lik(traits.a.2$values+traits.a.2$errors,
        params.a.2) -
      patherit::lik.poumm(traits.a.2$values[,1]+traits.a.2$errors[,1],
                          params.a.2$tree,
                          params.a.2$A[1,1,1],
                          params.a.2$Theta[1,1],
                          sqrt(params.a.2$Sigma[1,1,1]),
                          sqrt(params.a.2$Sigmae[1,1,1]),
                          distgr='maxlik', usempfr=0)) < EPS)
  expect_true(abs(
    lik(traits.a.3$values+traits.a.3$errors,
        params.a.3) -
      patherit::lik.poumm(traits.a.3$values[,1]+traits.a.3$errors[,1],
                          params.a.3$tree,
                          params.a.3$A[1,1,1],
                          params.a.3$Theta[1,1],
                          sqrt(params.a.3$Sigma[1,1,1]),
                          sqrt(params.a.3$Sigmae[1,1,1]),
                          distgr='maxlik', usempfr=0)) < EPS)
})

context(ctx <- "R=1/k=1/N=400")

# number of tips
N <- 400

# tree with one regime
tree.a <- phytools::pbtree(n=N, scale=1)

params.a.1$tree <- params.a.2$tree <- params.a.3$tree <- tree.a

# generate traits
traits.a.1 <- generateMVTrait(params.a.1, 0, verbose=TRUE)
traits.a.2 <- generateMVTrait(params.a.2, 0, verbose=TRUE)
traits.a.3 <- generateMVTrait(params.a.2, 0, verbose=TRUE)

test_that(paste(ctx, "Match univariate likelihood from patherit regime a"), {
  expect_true(
    abs(lik(traits.a.1$values+traits.a.1$errors,
            params.a.1) -
          patherit::lik.poumm(traits.a.1$values[,1]+traits.a.1$errors[,1],
                              params.a.1$tree,
                              params.a.1$A[1,1,1],
                              params.a.1$Theta[1,1],
                              sqrt(params.a.1$Sigma[1,1,1]),
                              sqrt(params.a.1$Sigmae[1,1,1]),
                              distgr='maxlik', usempfr=0)) < EPS)
  expect_true(abs(
    lik(traits.a.2$values+traits.a.2$errors,
        params.a.2) -
      patherit::lik.poumm(traits.a.2$values[,1]+traits.a.2$errors[,1],
                          params.a.2$tree,
                          params.a.2$A[1,1,1],
                          params.a.2$Theta[1,1],
                          sqrt(params.a.2$Sigma[1,1,1]),
                          sqrt(params.a.2$Sigmae[1,1,1]),
                          distgr='maxlik', usempfr=0)) < EPS)
  expect_true(abs(
    lik(traits.a.3$values+traits.a.3$errors,
        params.a.3) -
      patherit::lik.poumm(traits.a.3$values[,1]+traits.a.3$errors[,1],
                          params.a.3$tree,
                          params.a.3$A[1,1,1],
                          params.a.3$Theta[1,1],
                          sqrt(params.a.3$Sigma[1,1,1]),
                          sqrt(params.a.3$Sigmae[1,1,1]),
                          distgr='maxlik', usempfr=0)) < EPS)
})


context(ctx <- "R=1/k=3/N=2")

params.a.123$tree <- tree.a

traits.a.123 <- generateMVTrait(params.a.123, c(0,0,0), verbose=TRUE)

test_that(paste(ctx, "Match multivariate likelihood of independent traits regime a"), {
  expect_true(
    abs(lik(traits.a.123$values+traits.a.123$errors, params.a.123) -
          (patherit::lik.poumm(traits.a.123$values[,1]+traits.a.123$errors[,1],
                              params.a.123$tree,
                              params.a.123$A[1,1,1],
                              params.a.123$Theta[1,1],
                              sqrt(params.a.123$Sigma[1,1,1]),
                              sqrt(params.a.123$Sigmae[1,1,1]),
                              distgr='maxlik', usempfr=0)+
          patherit::lik.poumm(traits.a.123$values[,2]+traits.a.123$errors[,2],
                              params.a.123$tree,
                              params.a.123$A[1,2,2],
                              params.a.123$Theta[1,2],
                              sqrt(params.a.123$Sigma[1,2,2]),
                              sqrt(params.a.123$Sigmae[1,2,2]),
                              distgr='maxlik', usempfr=0)+
          patherit::lik.poumm(traits.a.123$values[,3]+traits.a.123$errors[,3],
                              params.a.123$tree,
                              params.a.123$A[1,3,3],
                              params.a.123$Theta[1,3],
                              sqrt(params.a.123$Sigma[1,3,3]),
                              sqrt(params.a.123$Sigmae[1,3,3]),
                              distgr='maxlik', usempfr=0))) < EPS)
})



context(ctx <- "R=1/k=3/N=400")

params.a.123$tree <- tree.a

traits.a.123 <- generateMVTrait(params.a.123, c(0,0,0), verbose=TRUE)


## Calculate likelihood
test_that(paste(ctx, "Match multivariate likelihood of independent traits regime a"), {
  expect_true(
    abs(lik(traits.a.123$values+traits.a.123$errors, params.a.123) -
          (patherit::lik.poumm(traits.a.123$values[,1]+traits.a.123$errors[,1],
                               params.a.123$tree,
                               params.a.123$A[1,1,1],
                               params.a.123$Theta[1,1],
                               sqrt(params.a.123$Sigma[1,1,1]),
                               sqrt(params.a.123$Sigmae[1,1,1]),
                               distgr='maxlik', usempfr=0)+
             patherit::lik.poumm(traits.a.123$values[,2]+traits.a.123$errors[,2],
                                 params.a.123$tree,
                                 params.a.123$A[1,2,2],
                                 params.a.123$Theta[1,2],
                                 sqrt(params.a.123$Sigma[1,2,2]),
                                 sqrt(params.a.123$Sigmae[1,2,2]),
                                 distgr='maxlik', usempfr=0)+
             patherit::lik.poumm(traits.a.123$values[,3]+traits.a.123$errors[,3],
                                 params.a.123$tree,
                                 params.a.123$A[1,3,3],
                                 params.a.123$Theta[1,3],
                                 sqrt(params.a.123$Sigma[1,3,3]),
                                 sqrt(params.a.123$Sigmae[1,3,3]),
                                 distgr='maxlik', usempfr=0))) < EPS)
})


plot(traits.a.1$values)


params.b.123$tree <- tree.a

traits.b.123 <- generateMVTrait(params.b.123, c(0,0,0), verbose=TRUE)

tree.ab <- phytools::sim.history(tree.a, Q, anc='a')

# convert the simmap tree to a normal phylo object with singleton nodes at the
# within-branch regime changes. The regimes are encoded as names of the edge.length
# vector
tree.ab.singles <- map.to.singleton(tree.ab)

params.ab.123$tree <- tree.ab



traits.ab.123 <- generateMVTrait(params.ab.123, c(0,0,0), verbose=TRUE)


