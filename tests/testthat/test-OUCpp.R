library(phytools)
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

## Specifying a bivariate OU process for each regime
# First, specify the A, theta, Sigma and sigmae2 parameters for each regime.
# Then we use the abind function to stack the parameters into arrays which's first
# dimension is the regime

# regimes

# in regime 'a' the three traits evolve according to three independent OU processes
a.X0 <- c(5, 2, 1)
a.Alpha <- rbind(
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
b.Alpha <- rbind(
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

Alpha <- abind::abind(a.Alpha, b.Alpha, along=-1, new.names=list(regime=c('a','b'), x=NULL, y=NULL))
Theta <- abind::abind(a.Theta, b.Theta, along=-1, new.names=list(regime=c('a', 'b'), xy=NULL))
Sigma <- abind::abind(a.Sigma, b.Sigma, along=-1, new.names=list(regime=c('a','b'), x=NULL, y=NULL))
Sigmae <- abind::abind(a.Sigmae2, b.Sigmae2, along=-1, new.names=list(regime=c('a','b'), x=NULL, y=NULL))
# regimes 'a' and 'b', traits 1, 2 and 3
model.ab.123 <- list(X0 = a.X0,
                     Alpha=Alpha[,,,drop=FALSE],
                     Theta=Theta[,,drop=FALSE],
                     Sigma=Sigma[,,,drop=FALSE],
                     Sigmae=Sigmae[,,,drop=FALSE])
class(model.ab.123) <- 'OU'

tree.a <- rtree(3)

tree.ab <- phytools::sim.history(tree.a, Q, anc='a')

# convert the simmap tree to a normal phylo object with singleton nodes at the
# within-branch regime changes. The regimes are encoded as names of the edge.length
# vector
tree.ab.singles <- map.to.singleton(tree.ab)

tree.ab.singles$edge.regime <- names(tree.ab.singles$edge.length)

traits.ab.123 <- mvsim(tree.ab.singles, model.ab.123, c(0,0,0), verbose=TRUE)

mvlik(traits.ab.123$values+traits.ab.123$errors, tree.ab.singles, model.ab.123)


X <- t(traits.ab.123$values)
Pc <- apply(X, 1:2, function(x) as.integer(!is.na(x)))

objectOU <- PCMBase:::QuadraticPolynomialOU$new(tree = tree.ab.singles,
                                                regimes_unique = c("a","b"),
                                                X = X[,1:length(tree.ab.singles$tip.label)],
                                                Pc = Pc[,1:length(tree.ab.singles$tip.label)])

llR <- mvlik(X = traits.ab.123$values, tree = tree.ab.singles, model = model.ab.123, verbose = TRUE)

llCpp <- mvlik(tree=tree.ab.singles, model = model.ab.123, pruneI = objectOU, pc = NULL)
