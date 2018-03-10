library(ape)
library(testthat)
library(mvtnorm)
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

# regimes

# in regime 'a' the three traits evolve according to three independent OU processes
a.X0 <- c(5, 2, 1)
a.H <- rbind(c(0, 0, 0),
                 c(0, 2, 0),
                 c(0, 0, 3))
a.Theta <- c(10, 6, 2)
a.Sigma <- rbind(c(1.6, 0, 0),
                 c(0, 2.4, 0),
                 c(0, 0, 2))
a.Sigmae2 <- rbind(c(0, 0, 0),
                   c(0, 0, 0),
                   c(0, 0, 0))
a.Sigmaj <- rbind(c(.2, 0, 0),
                  c(0, .3, 0),
                  c(0, 0, .4))
a.mj <- c(1, 4, 7)

# in regime 'b' there is correlation between the traits
b.X0 <- c(12, 4, 3)
b.H <- rbind(c(2, .1, .2),
                 c(.1, .6, .2),
                 c(.2, .2, .3))
b.Theta <- c(10, 6, 2)
b.Sigma <- rbind(c(1.6, .3, .3),
                 c(.3, 0.3, .4),
                 c(.3, .4, 2))
b.Sigmae2 <- rbind(c(.2, 0, 0),
                   c(0, .3, 0),
                   c(0, 0, .4))
b.Sigmaj <- rbind(c(.2, 0.1, 0.1),
                  c(0.1, .3, 0),
                  c(0.1, 0, .4))
b.mj <- c(11, 17, 42)

# First, specify the H, theta, Sigma sigmae2, Sigmaj, mj parameters for each regime.
# Then we use the abind function to stack the parameters into arrays which's first
# dimension is the regime

H <- abind(a.H, b.H, along=3, new.names=list(x=NULL, y=NULL, regime=c('a','b')))
Theta <- abind(a.Theta, b.Theta, along=2, new.names=list(xy=NULL, regime=c('a', 'b')))
Sigma <- abind(a.Sigma, b.Sigma, along=3, new.names=list(x=NULL, y=NULL, regime=c('a','b')))
Sigmae <- abind(a.Sigmae2, b.Sigmae2, along=3, new.names=list(x=NULL, y=NULL, regime=c('a','b')))
Sigmaj <- abind(a.Sigmaj, b.Sigmaj, along=3, new.names=list(x=NULL, y=NULL, regime=c('a','b')))
mj <- abind(a.mj, b.mj, along=2, new.names=list(xy=NULL, regime=c('a', 'b')))


# regime 'a', traits 1, 2 and 3
model.a.123 <- PCM("JOU", 3, "a", list(X0 = a.X0,
                                       H=H[,,'a',drop=FALSE],
                                       Theta=Theta[,'a',drop=FALSE],
                                       Sigma=Sigma[,,'a',drop=FALSE],
                                       Sigmae=Sigmae[,,'a',drop=FALSE],
                                       Sigmaj=Sigmaj[,,'a',drop=FALSE],
                                       mj=mj[,'a',drop=FALSE]))

# regime 'b', traits 1, 2 and 3
model.b.123 <- PCM("JOU", 3, "b", list(X0 = b.X0,
                                       H=H[,,'b',drop=FALSE],
                                       Theta=Theta[,'b',drop=FALSE],
                                       Sigma=Sigma[,,'b',drop=FALSE],
                                       Sigmae=Sigmae[,,'b',drop=FALSE],
                                       Sigmaj=Sigmaj[,,'b',drop=FALSE],
                                       mj=mj[,'b',drop=FALSE]))

# regimes 'a' and 'b', traits 1, 2 and 3
model.ab.123 <- PCM("JOU", 3, c("a", "b"), list(X0 = a.X0,
                                                H=H[,,,drop=FALSE],
                                                Theta=Theta[,,drop=FALSE],
                                                Sigma=Sigma[,,,drop=FALSE],
                                                Sigmae=Sigmae[,,,drop=FALSE],
                                                Sigmaj=Sigmaj[,,,drop=FALSE],
                                                mj=mj[,,drop=FALSE]))


context(ctx <- "R=1/k=1/N=5")

# number of tips
N <- 400

# tree with one regime

tree.a <- rtree(N) #phytools::pbtree(n=N, scale=1)
tree.b <- rtree(N) #phytools::pbtree(n=N, scale=1)


tree.a$edge.jump <- sample(as.integer(0:1), size = nrow(tree.a$edge), replace = TRUE)
tree.b$edge.jump <- tree.a$edge.jump

if(require(phytools)) {
  # tree with two regimes
  tree.ab <- phytools::sim.history(tree.a, Q, anc='a')

  # convert the simmap tree to a normal phylo object with singleton nodes at the
  # within-branch regime changes. The regimes are encoded as names of the edge.length
  # vector
  tree.ab.singles <- map.to.singleton(tree.ab)
  tree.ab.singles$edge.regime <- names(tree.ab.singles$edge.length)
  # generate traits
} else {
  tree.ab <- tree.a
  tree.ab$edge.regime <- sample(c("a", "b"), size = length(tree.ab$edge.length), replace = TRUE)
  tree.ab.singles <- tree.ab
}

tree.ab.singles$edge.jump = sample(as.integer(0:1), size = nrow(tree.ab.singles$edge), replace = TRUE)

traits.a.123 <- PCMSim(tree.a, model.a.123, c(0,0,0), verbose=TRUE)
traits.b.123 <- PCMSim(tree.b, model.b.123, c(0,0,0), verbose=TRUE)

# The singles tree is used for the generation of the traits in this case. Since the singles tree will be used
# for the likelihood calculation, the length of xi has to match the number of edges in the singles tree. Therefore
# in case the tree.ab is used for the generation of the traits then the length of xi will not match the number of
# edges which is less in this tree.
traits.ab.123 <- PCMSim(tree.ab.singles, model.ab.123, c(0,0,0), verbose=TRUE)


# calculate likelihoods

JOU.lik.a <-  PCMLik(traits.a.123$values+traits.a.123$errors, tree.a, model.a.123)
JOU.lik.b <-  PCMLik(traits.b.123$values+traits.b.123$errors, tree.b, model.b.123)
JOU.lik.ab <-  PCMLik(traits.ab.123$values+traits.ab.123$errors, tree.ab.singles, model.ab.123)



if(require(PCMBaseCpp)) {
  cat("Testing PCMBaseCpp on JOU:\n")

  test_that("a.123",
            expect_equal(PCMLik(traits.a.123$values+traits.a.123$errors, tree.a, model.a.123),
                         PCMLik(traits.a.123$values+traits.a.123$errors, tree.a, model.a.123,
                               pruneI = PCMCppPruningObject(X = traits.a.123$values[, 1:length(tree.a$tip.label)],
                                                     tree = tree.a,
                                                     model.a.123))))


  cat("Testing PCMBaseCpp on JOU with missing values:\n")

  values <- traits.ab.123$values[, 1:length(tree.ab.singles$tip.label)] +
    traits.ab.123$errors[, 1:length(tree.ab.singles$tip.label)]

  pruneI <- PCMCppPruningObject(X = values, tree = tree.ab.singles, model.ab.123)

  PCMAbCdEf(tree.ab.singles, model.ab.123,
            PCMValidate(tree.ab.singles, model.ab.123),
            PCMPresentCoordinates(values, tree.ab.singles))

  test_that("ab.123",
            expect_equal(PCMLik(values, tree.ab.singles, model.ab.123),
                         PCMLik(tree = tree.ab.singles, model = model.ab.123, pruneI = pruneI)))

}

