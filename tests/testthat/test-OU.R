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
model.a.1 <- list(X0 = a.X0[1],
                  H=H[1,1,'a',drop=FALSE],
                  Theta=Theta[1,'a',drop=FALSE],
                  Sigma=Sigma[1,1,'a',drop=FALSE],
                  Sigmae=Sigmae[1,1,'a',drop=FALSE])
class(model.a.1) <- 'OU'

# regime 'a', trait 2
model.a.2 <- list(X0 = a.X0[2],
                  H=H[2,2,'a',drop=FALSE],
                  Theta=Theta[2,'a',drop=FALSE],
                  Sigma=Sigma[2,2,'a',drop=FALSE],
                  Sigmae=Sigmae[2,2,'a',drop=FALSE])
class(model.a.2) <- 'OU'

# regime 'a', trait 3
model.a.3 <- list(X0 = a.X0[3],
                  H=H[3,3,'a',drop=FALSE],
                  Theta=Theta[3,'a',drop=FALSE],
                  Sigma=Sigma[3,3,'a',drop=FALSE],
                  Sigmae=Sigmae[3,3,'a',drop=FALSE])
class(model.a.3) <- 'OU'

# regime 'a', traits 1, 2 and 3
model.a.123 <- list(X0 = a.X0,
                    H=H[,,'a',drop=FALSE],
                    Theta=Theta[,'a',drop=FALSE],
                    Sigma=Sigma[,,'a',drop=FALSE],
                    Sigmae=Sigmae[,,'a',drop=FALSE])
class(model.a.123) <- 'OU'


# regime 'b', traits 1, 2 and 3
model.b.123 <- list(X0 = b.X0,
                    H=H[,,'b',drop=FALSE],
                    Theta=Theta[,'b',drop=FALSE],
                    Sigma=Sigma[,,'b',drop=FALSE],
                    Sigmae=Sigmae[,,'b',drop=FALSE])
class(model.b.123) <- 'OU'

# regimes 'a' and 'b', traits 1, 2 and 3
model.ab.123 <- list(X0 = a.X0,
                     H=H[,,,drop=FALSE],
                     Theta=Theta[,,drop=FALSE],
                     Sigma=Sigma[,,,drop=FALSE],
                     Sigmae=Sigmae[,,,drop=FALSE])
class(model.ab.123) <- 'OU'


context(ctx <- "R=1/k=1/N=2")

# number of tips
N <- 2

# tree with one regime
tree.a <- rtree(N) # phytools::pbtree(n=N, scale=1)
tree.a$edge.regime <- rep("a", length(tree.a$edge.length))

# generate traits
traits.a.1 <- PCMSim(tree.a, model.a.1, 0, verbose=TRUE)
traits.a.2 <- PCMSim(tree.a, model.a.2, 0, verbose=TRUE)
traits.a.3 <- PCMSim(tree.a, model.a.2, 0, verbose=TRUE)


# test likelihood
test_that(paste(ctx, "Match univariate likelihood from patherit regime a"), {
  expect_true(
    abs(PCMLik(traits.a.1$values + traits.a.1$errors, tree.a,
            model.a.1) -
          POUMM::likPOUMMGivenTreeVTips(
            z = traits.a.1$values[1,] + traits.a.1$errors[1,],
            tree = tree.a,
            alpha = model.a.1$H[1,1,1],
            theta = model.a.1$Theta[1,1],
            sigma = sqrt(model.a.1$Sigma[1,1,1]),
            sigmae = sqrt(model.a.1$Sigmae[1,1,1]),
            g0 = a.X0[1])
        ) < EPS)
  expect_true(
    abs(PCMLik(traits.a.2$values + traits.a.2$errors, tree.a,
        model.a.2) -
          POUMM::likPOUMMGivenTreeVTips(
            z = traits.a.2$values[1,] + traits.a.2$errors[1,],
            tree = tree.a,
            alpha = model.a.2$H[1,1,1],
            theta = model.a.2$Theta[1,1],
            sigma = sqrt(model.a.2$Sigma[1,1,1]),
            sigmae = sqrt(model.a.2$Sigmae[1,1,1]),
            g0 = a.X0[2])
        ) < EPS)
  expect_true(
    abs(PCMLik(traits.a.3$values + traits.a.3$errors, tree.a,
        model.a.3) -
          POUMM::likPOUMMGivenTreeVTips(
            z = traits.a.3$values[1,] + traits.a.3$errors[1,],
            tree = tree.a,
            alpha = model.a.3$H[1,1,1],
            theta = model.a.3$Theta[1,1],
            sigma = sqrt(model.a.3$Sigma[1,1,1]),
            sigmae = sqrt(model.a.3$Sigmae[1,1,1]),
            g0 = a.X0[3])) < EPS)
})

context(ctx <- "R=1/k=1/N=400")

# number of tips
N <- 400

# tree with one regime
tree.a <- rtree(N) # pbtree(n=N, scale=1)
tree.a$edge.regime <- rep("a", length(tree.a$edge.length))

# generate traits
traits.a.1 <- PCMSim(tree.a, model.a.1, 0, verbose=TRUE)
traits.a.2 <- PCMSim(tree.a, model.a.2, 0, verbose=TRUE)
traits.a.3 <- PCMSim(tree.a, model.a.2, 0, verbose=TRUE)

test_that(paste(ctx, "Match univariate likelihood from patherit regime a"), {
  expect_true(
    abs(PCMLik(traits.a.1$values + traits.a.1$errors, tree.a,
              model.a.1) -
          POUMM::likPOUMMGivenTreeVTips(
            z = traits.a.1$values[1,] + traits.a.1$errors[1,],
            tree = tree.a,
            alpha = model.a.1$H[1,1,1],
            theta = model.a.1$Theta[1,1],
            sigma = sqrt(model.a.1$Sigma[1,1,1]),
            sigmae = sqrt(model.a.1$Sigmae[1,1,1]),
            g0 = a.X0[1])
    ) < EPS)
  expect_true(
    abs(PCMLik(traits.a.2$values + traits.a.2$errors, tree.a,
              model.a.2) -
          POUMM::likPOUMMGivenTreeVTips(
            z = traits.a.2$values[1,] + traits.a.2$errors[1,],
            tree = tree.a,
            alpha = model.a.2$H[1,1,1],
            theta = model.a.2$Theta[1,1],
            sigma = sqrt(model.a.2$Sigma[1,1,1]),
            sigmae = sqrt(model.a.2$Sigmae[1,1,1]),
            g0 = a.X0[2])
    ) < EPS)
  expect_true(
    abs(PCMLik(traits.a.3$values + traits.a.3$errors, tree.a,
              model.a.3) -
          POUMM::likPOUMMGivenTreeVTips(
            z = traits.a.3$values[1,] + traits.a.3$errors[1,],
            tree = tree.a,
            alpha = model.a.3$H[1,1,1],
            theta = model.a.3$Theta[1,1],
            sigma = sqrt(model.a.3$Sigma[1,1,1]),
            sigmae = sqrt(model.a.3$Sigmae[1,1,1]),
            g0 = a.X0[3])) < EPS)
})

context(ctx <- "R=1/k=3/N=2")

traits.a.123 <- PCMSim(tree.a, model.a.123, c(0,0,0), verbose=TRUE)

test_that(paste(ctx, "Match multivariate likelihood of independent traits regime a"), {
  expect_true(
    abs(PCMLik(traits.a.123$values+traits.a.123$errors, tree.a, model.a.123) -
          (POUMM::likPOUMMGivenTreeVTips(
            traits.a.123$values[1,]+traits.a.123$errors[1,],
            tree.a,
            model.a.123$H[1,1,1],
            model.a.123$Theta[1,1],
            sqrt(model.a.123$Sigma[1,1,1]),
            sqrt(model.a.123$Sigmae[1,1,1]),
            a.X0[1]) +

             POUMM::likPOUMMGivenTreeVTips(
               traits.a.123$values[2,]+traits.a.123$errors[2,],
               tree.a,
               model.a.123$H[2,2,1],
               model.a.123$Theta[2,1],
               sqrt(model.a.123$Sigma[2,2,1]),
               sqrt(model.a.123$Sigmae[2,2,1]),
               a.X0[2]) +

             POUMM::likPOUMMGivenTreeVTips(traits.a.123$values[3,]+traits.a.123$errors[3,],
                              tree.a,
                              model.a.123$H[3,3,1],
                              model.a.123$Theta[3,1],
                              sqrt(model.a.123$Sigma[3,3,1]),
                              sqrt(model.a.123$Sigmae[3,3,1]),
                              a.X0[3]))
        ) < EPS)
})



context(ctx <- "R=1/k=3/N=400")

traits.a.123 <- PCMSim(tree.a, model.a.123, c(0,0,0), verbose=TRUE)


## Calculate likelihood
test_that(paste(ctx, "Match multivariate likelihood of independent traits regime a"), {
  expect_true(
    abs(PCMLik(traits.a.123$values+traits.a.123$errors, tree.a, model.a.123) -
          (POUMM::likPOUMMGivenTreeVTips(
            traits.a.123$values[1,]+traits.a.123$errors[1,],
            tree.a,
            model.a.123$H[1,1,1],
            model.a.123$Theta[1,1],
            sqrt(model.a.123$Sigma[1,1,1]),
            sqrt(model.a.123$Sigmae[1,1,1]),
            a.X0[1]) +
             POUMM::likPOUMMGivenTreeVTips(
               traits.a.123$values[2,] + traits.a.123$errors[2,],
               tree.a,
               model.a.123$H[2,2,1],
               model.a.123$Theta[2,1],
               sqrt(model.a.123$Sigma[2,2,1]),
               sqrt(model.a.123$Sigmae[2,2,1]),
               a.X0[2]) +
             POUMM::likPOUMMGivenTreeVTips(
               traits.a.123$values[3,]+traits.a.123$errors[3,],
               tree.a,
               model.a.123$H[3,3,1],
               model.a.123$Theta[3,1],
               sqrt(model.a.123$Sigma[3,3,1]),
               sqrt(model.a.123$Sigmae[3,3,1]),
               a.X0[3]))
        ) < EPS)
})

plot(traits.a.1$values)


traits.b.123 <- PCMSim(tree.a, model.b.123, c(0,0,0), verbose=TRUE)

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

traits.ab.123 <- PCMSim(tree.ab.singles, model.ab.123, c(0,0,0), verbose=TRUE)

if(require(PCMBaseCpp)) {
  cat("Testing PCMBaseCpp on OU:\n")

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
  values[sample(x=1:length(values), 88)] <- NA

  pruneInfoR <- PCMPruningOrder(tree.ab.singles)
  pruneInfoCpp <- PCMCppPruningObject(X = values,
                               tree = tree.ab.singles,
                               model.ab.123)
  test_that("ab.123 with missing values",
            expect_equal(PCMLik(values, tree.ab.singles, model.ab.123,
                               pruneI = pruneInfoR),
                         PCMLik(values, tree.ab.singles, model.ab.123,
                               pruneI = pruneInfoCpp)))

  print(PCMLik(traits.ab.123$values + traits.ab.123$errors, tree.ab.singles, model.ab.123,
        pruneI = PCMCppPruningObject(X = traits.ab.123$values[, 1:length(tree.ab.singles$tip.label)] +
                                traits.ab.123$errors[, 1:length(tree.ab.singles$tip.label)],
                              tree = tree.ab.singles,
                              model.ab.123)))

  print(PCMLik(values, tree.ab.singles, model.ab.123,
              pruneI = PCMCppPruningObject(X = values,
                                    tree = tree.ab.singles,
                                    model.ab.123)))


  if(require(microbenchmark)) {
    cat("microbenchmark test")

    options(PCMBase.PCMLmr.mode=11)
    print(microbenchmark(
      PCMLik(values, tree.ab.singles, model.ab.123, pruneI = pruneInfoR),
      PCMLik(values, tree.ab.singles, model.ab.123, pruneI = pruneInfoCpp),
      times = 10
    ))

    options(PCMBase.PCMLmr.mode=21)
    print(microbenchmark(
      PCMLik(values, tree.ab.singles, model.ab.123, pruneI = pruneInfoCpp)
    ))
  }
}

