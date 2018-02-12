library(ape)
library(testthat)
library(PCMBase)


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

# regimes

# in regime 'a' the three traits evolve according to three independent OU processes
a.X0 <- c(5, 2, 1)
a.H <- rbind(c(0, 0, 0),
                 c(0, 2, 0),
                 c(0, 0, 3))
a.H2 <- rbind(c(0, 0, 0),
                 c(0, 2, 0),
                 c(0, 0, 3))

a.Theta <- c(10, 6, 2)
a.Sigma <- rbind(c(1.6, 0, 0),
                 c(0, 2.4, 0),
                 c(0, 0, 2))
a.Sigmae2 <- rbind(c(0, 0, 0),
                   c(0, 0, 0),
                   c(0, 0, 0))

# in regime 'b' there is correlation between the traits
b.X0 <- c(12, 4, 3)
b.H <- rbind(c(2, .1, .2),
                 c(.1, .6, .2),
                 c(.2, .2, .3))

b.H2 <- rbind(c(2, .1, .2),
                 c(.1, .6, .2),
                 c(.2, .2, .3))

b.Theta <- c(10, 6, 2)
b.Sigma <- rbind(c(1.6, .3, .3),
                 c(.3, 0.3, .4),
                 c(.3, .4, 2))
b.Sigmae2 <- rbind(c(.2, 0, 0),
                   c(0, .3, 0),
                   c(0, 0, .4))

# First, specify the ALpha1,H2, theta, Sigma and sigmae2 parameters for each regime.
# Then we use the abind function to stack the parameters into arrays which's first
# dimension is the regime


H1 <- abind::abind(a.H, b.H, along=3, new.names=list(x=NULL, y=NULL, regime=c('a','b')))
H2 <- abind::abind(a.H2, b.H2, along=3, new.names=list(x=NULL, y=NULL, regime=c('a','b')))
Theta <- abind::abind(a.Theta, b.Theta, along=2, new.names=list(xy=NULL, regime=c('a', 'b')))
Sigma <- abind::abind(a.Sigma, b.Sigma, along=3, new.names=list(x=NULL, y=NULL, regime=c('a','b')))
Sigmae <- abind::abind(a.Sigmae2, b.Sigmae2, along=3, new.names=list(x=NULL, y=NULL, regime=c('a','b')))


# regime 'a', traits 1, 2 and 3
model.a.123 <- list(X0 = a.X0,
                    H=H1[,,'a',drop=FALSE],
                    Theta=Theta[,'a',drop=FALSE],
                    Sigma=Sigma[,,'a',drop=FALSE],
                    Sigmae=Sigmae[,,'a',drop=FALSE])
class(model.a.123) <- 'OU'

model.a.123.TwoSpeedOU <- list(X0 = a.X0,
                    H1=H1[,,'a',drop=FALSE],
                    H2=H2[,,'a',drop=FALSE],
                    Theta=Theta[,'a',drop=FALSE],
                    Sigma=Sigma[,,'a',drop=FALSE],
                    Sigmae=Sigmae[,,'a',drop=FALSE])
class(model.a.123.TwoSpeedOU) <- 'TwoSpeedOU'

#####################################################################################################

context(ctx <- "R=1/k=3/N=400")

# number of tips
N <- 400

# tree with one regime
tree.a <- rtree(N) # phytools::pbtree(n=N, scale=1)
tree.a$edge.regime <- names(tree.a$edge.length)

# generate traits

traits.a.123 <- PCMSim(tree.a, model.a.123, c(0,0,0), verbose=TRUE)

## Calculate likelihood

lik.OU = PCMLik(traits.a.123$values+traits.a.123$errors, tree.a, model.a.123)

lik.TwoSpeedOU = PCMLik(traits.a.123$values+traits.a.123$errors, tree.a, model.a.123.TwoSpeedOU)

cat('OU likelihood=',lik.OU,'\n')
cat('TwoSpeedOU likelihood=',lik.TwoSpeedOU,'\n')

test_that(paste(ctx, "Match multivariate likelihood of independent traits regime a"), {
  expect_true(abs(lik.OU - lik.TwoSpeedOU) < EPS)
})

#####################################################################################################

# regime 'a', traits 1, 2 and 3
model.b.123 <- list(X0 = b.X0,
                    H=H1[,,'b',drop=FALSE],
                    Theta=Theta[,'b',drop=FALSE],
                    Sigma=Sigma[,,'b',drop=FALSE],
                    Sigmae=Sigmae[,,'b',drop=FALSE])
class(model.b.123) <- 'OU'

model.b.123.TwoSpeedOU <- list(X0 = b.X0,
                          H1=H1[,,'b',drop=FALSE],
                          H2=H2[,,'b',drop=FALSE],
                          Theta=Theta[,'b',drop=FALSE],
                          Sigma=Sigma[,,'b',drop=FALSE],
                          Sigmae=Sigmae[,,'b',drop=FALSE])
class(model.b.123.TwoSpeedOU) <- 'TwoSpeedOU'

context(ctx <- "R=1/k=3/N=400")

# number of tips
N <- 400

# tree with one regime
tree.b <- rtree(N) # phytools::pbtree(n=N, scale=1)
tree.b$edge.regime <- names(tree.b$edge.length)

# generate traits

traits.b.123 <- PCMSim(tree.b, model.b.123, c(0,0,0), verbose=TRUE)

## Calculate likelihood

lik.OU = PCMLik(traits.b.123$values+traits.b.123$errors, tree.b, model.b.123)

lik.TwoSpeedOU = PCMLik(traits.b.123$values+traits.b.123$errors, tree.b, model.b.123.TwoSpeedOU)

cat('OU likelihood=',lik.OU,'\n')
cat('TwoSpeedOU likelihood=',lik.TwoSpeedOU,'\n')

test_that(paste(ctx, "Match multivariate likelihood of dependent traits regime b"), {
  expect_true(abs(lik.OU - lik.TwoSpeedOU) < EPS)
})


if(require(phytools)) {
  tree.ab <- phytools::sim.history(tree.a, Q, anc='a')
  tree.ab$edge.regime <- names(tree.ab$edge.length)

  # convert the simmap tree to a normal phylo object with singleton nodes at the
  # within-branch regime changes. The regimes are encoded as names of the edge.length
  # vector
  tree.ab.singles <- map.to.singleton(tree.ab)
  tree.ab.singles$edge.regime <- names(tree.ab.singles$edge.length)
} else {
  tree.ab <- tree.a
  tree.ab$edge.regime <- sample(c("a", "b"), size = length(tree.ab$edge.length), replace = TRUE)
  tree.ab.singles <- tree.ab
}

model.ab.123 <- list(X0 = a.X0,
                     H1=H1[,,,drop=FALSE],
                     H2=H2[,,,drop=FALSE],
                     Theta=Theta[,,drop=FALSE],
                     Sigma=Sigma[,,,drop=FALSE],
                     Sigmae=Sigmae[,,,drop=FALSE])
class(model.ab.123) <- 'TwoSpeedOU'

traits.ab.123 <- PCMSim(tree.ab.singles, model.ab.123, c(0,0,0), verbose=TRUE)


if(require(PCMBaseCpp)) {
  cat("Testing PCMBaseCpp on TwoSpeedOU:\n")

  test_that("a.123",
            expect_equal(PCMLik(traits.a.123$values+traits.a.123$errors, tree.a, model.a.123),
                         PCMLik(traits.a.123$values+traits.a.123$errors, tree.a, model.a.123,
                               pruneI = PCMCppPruningObject(X = traits.a.123$values[, 1:length(tree.a$tip.label)],
                                                     tree = tree.a,
                                                     model.a.123))))

  test_that("ab.123",
            expect_equal(PCMLik(traits.ab.123$values + traits.ab.123$errors, tree.ab.singles, model.ab.123),
                         PCMLik(traits.ab.123$values + traits.ab.123$errors, tree.ab.singles, model.ab.123,
                               pruneI = PCMCppPruningObject(X = traits.ab.123$values[, 1:length(tree.ab.singles$tip.label)] +
                                                       traits.ab.123$errors[, 1:length(tree.ab.singles$tip.label)],
                                                     tree = tree.ab.singles,
                                                     model.ab.123))))


  values <- traits.ab.123$values[, 1:length(tree.ab.singles$tip.label)] + traits.ab.123$errors[, 1:length(tree.ab.singles$tip.label)]
  values[sample(x=1:length(values), 50)] <- NA

  pruneInfoR <- PCMPruningOrder(tree.ab.singles)
  pruneInfoCpp <- PCMCppPruningObject(X = values,
                               tree = tree.ab.singles,
                               model.ab.123)

  test_that("ab.123 with missing values",
            expect_equal(PCMLik(values, tree.ab.singles, model.ab.123, pruneI = pruneInfoR),
                         PCMLik(values, tree.ab.singles, model.ab.123,
                               pruneI = pruneInfoCpp)))


  print(PCMLik(values, tree.ab.singles, model.ab.123, pruneI = pruneInfoCpp))

  print(PCMLik(values, tree.ab.singles, model.ab.123,
              pruneI = PCMCppPruningObject(X = values,
                                    tree = tree.ab.singles,
                                    model.ab.123)))


  # if(require(microbenchmark)) {
  #   cat("microbenchmark test")
  #
  #   options(PCMBase.PCMLmr.mode=11)
  #   print(microbenchmark(
  #     PCMLik(values, tree.ab.singles, model.ab.123, pruneI = pruneInfoR),
  #     PCMLik(values, tree.ab.singles, model.ab.123, pruneI = pruneInfoCpp), times = 10
  #   ))
  #
  #   options(PCMBase.PCMLmr.mode=21)
  #   print(microbenchmark(
  #     PCMLik(values, tree.ab.singles, model.ab.123, pruneI = pruneInfoCpp)
  #   ))
  # }
}



