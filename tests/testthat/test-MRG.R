library(ape)
library(testthat)
library(PCMBase)
library(abind)

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

# in regime 'a' the three traits evolve according to three independent OU processes
a.X0 <- c(5, 2, 1)
a.H <- rbind(
  c(0, 0, 0),
  c(0, 2, 0),
  c(0, 0, 3))
a.Theta <- c(10, 6, 2)
a.Sigma <- rbind(
  c(1.6, 0, 0),
  c(0, 2.4, 0),
  c(0, 0, 2))
a.Sigmae2 <- rbind(
  c(0, 0, 0),
  c(0, 0, 0),
  c(0, 0, 0))

# in regime 'b' there is correlation between the traits
b.X0 <- c(12, 4, 3)
b.H <- rbind(
  c(2, .1, .2),
  c(.1, .6, .2),
  c(.2, .2, .3))
b.Theta <- c(10, 6, 2)
b.Sigma <- rbind(
  c(1.6, .3, .3),
  c(.3, 0.3, .4),
  c(.3, .4, 2))
b.Sigmae2 <- rbind(
  c(.2, 0, 0),
  c(0, .3, 0),
  c(0, 0, .4))

H <- abind(a.H, b.H, along=3, new.names=list(x=NULL, y=NULL, regime=c('a','b')))
Theta <- abind(a.Theta, b.Theta, along=2, new.names=list(xy=NULL, regime=c('a','b')))
Sigma <- abind(a.Sigma, b.Sigma, along=3, new.names=list(x=NULL, y=NULL, regime=c('a','b')))
Sigmae <- abind(a.Sigmae2, b.Sigmae2, along=3, new.names=list(x=NULL, y=NULL, regime=c('a','b')))


## Simulations of trait data

# regime 'a', trait 1
model.a.1 <- PCM("OU1", 1, "a", params = list(H=H[1,1,'a',drop=FALSE],
                                             Theta=Theta[1,'a',drop=FALSE],
                                             Sigma=Sigma[1,1,'a',drop=FALSE],
                                             Sigmae=Sigmae[1,1,'a',drop=FALSE]))


# regime 'a', traits 1, 2 and 3
model.a.123 <- PCM("OU1", 3, "a", params = list(H=H[,,'a',drop=FALSE],
                                               Theta=Theta[,'a',drop=FALSE],
                                               Sigma=Sigma[,,'a',drop=FALSE],
                                               Sigmae=Sigmae[,,'a',drop=FALSE]))


# regime 'b', traits 1, 2 and 3
model.b.123 <- PCM("OU1", 3, "b", params = list(H=H[,,'b',drop=FALSE],
                                               Theta=Theta[,'b',drop=FALSE],
                                               Sigma=Sigma[,,'b',drop=FALSE],
                                               Sigmae=Sigmae[,,'b',drop=FALSE]))

# regimes 'a' and 'b', traits 1, 2 and 3
model.ab.123 <- PCM("OU", 3, c("a", "b"), params = list(X0 = a.X0,
                                                        H=H[,,,drop=FALSE],
                                                        Theta=Theta[,,drop=FALSE],
                                                        Sigma=Sigma[,,,drop=FALSE],
                                                        Sigmae=Sigmae[,,,drop=FALSE]))

PCMParentClasses.MRG_ab <<- function(model) c("MRG", "GaussianPCM", "PCM")
PCMSpecifyParams.MRG_ab <<- function(model, ...) {
  k <- attr(model, "k")

  list(
    X0 = list(default = rep(0, k),
              type = c("gvector", "full"),
              description = "trait vector at the root; global for all model regimes"),
    'a' = list(
      default = PCM("OU1", k, 1),
      type = "model",
      description = "Stabilizing selection"),
    'b' = list(
      default = PCM("OU1", k, 1),
      type = c("model"),
      description = "Stabilizing selection"))
}

model_MRG_ab <- PCM("MRG_ab", k = 3, regimes = c("a", "b"))
PCMSetParams(model_MRG_ab, list(X0 = a.X0, a = model.a.123, b = model.b.123))


# number of tips
N <- 400

# tree with one regime
tree.a <- rtree(N) # pbtree(n=N, scale=1)
PCMSetDefaultRegime(tree.a, model.a.123)
#tree.a$edge.regime <- rep("a", length(tree.a$edge.length))

traits.a.123 <- PCMSim(tree.a, model.a.123, c(0,0,0), verbose=TRUE)


context(ctx <- "R=1/k=3/N=400")

traits.a.123 <- PCMSim(tree.a, model.a.123, c(0,0,0), verbose=TRUE)

PCMSetDefaultRegime(tree.a, model.b.123)
traits.b.123 <- PCMSim(tree.a, model.b.123, c(0,0,0), verbose=TRUE)

if(require(phytools)) {
  tree.ab <- phytools::sim.history(tree.a, Q, anc='a')
  tree.ab$edge.regime <- names(tree.ab$edge.length)


  # convert the simmap tree to a normal phylo object with singleton nodes at the
  # within-branch regime changes. The regimes are encoded as names of the edge.length
  # vector
  tree.ab.singles <- map.to.singleton(tree.ab)
  tree.ab.singles$edge.regime <- names(tree.ab.singles$edge.length)
  tree.ab.singles$edge.jump <- rep(0, length(tree.ab.singles$edge.length))
} else {
  tree.ab <- tree.a
  tree.ab$edge.regime <- sample(c("a", "b"), size = length(tree.ab$edge.length), replace = TRUE)
  tree.ab$edge.jump <- rep(0, length(tree.ab$edge.length))
  tree.ab.singles <- tree.ab
}

traits.ab.123 <- PCMSim(tree.ab.singles, model.ab.123, c(0,0,0), verbose=TRUE)
traits.ab.123 <- traits.ab.123[, 1:N]


test_that("Equal OU and MRG likelihoods", expect_equal(
  PCMLik(traits.ab.123, tree.ab.singles, model_MRG_ab),
  PCMLik(traits.ab.123, tree.ab.singles, model.ab.123)
))

if(require(PCMBaseCpp)) {
  cat("Testing PCMBaseCpp on MRG model:\n")
  test_that("Equal OU and MRG likelihoods", expect_equal(
    PCMLik(traits.ab.123, tree.ab.singles, model_MRG_ab,
           metaI = PCMInfoCpp(traits.ab.123, tree.ab.singles, model_MRG_ab)),
    PCMLik(traits.ab.123, tree.ab.singles, model.ab.123)
  ))

  metaI = PCMInfoCpp(traits.ab.123, tree.ab.singles, model_MRG_ab)
  microbenchmark::microbenchmark(PCMLik(traits.ab.123, tree.ab.singles, model_MRG_ab,
         metaI = metaI))
}

