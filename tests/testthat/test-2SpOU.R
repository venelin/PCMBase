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

# regimes

# in regime 'a' the three traits evolve according to three independent OU processes
a.X0 <- c(5, 2, 1)
a.Alpha <- rbind(c(0, 0, 0),
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
b.Alpha <- rbind(c(2, .1, .2),
                 c(.1, .6, .2),
                 c(.2, .2, .3))
b.Theta <- c(10, 6, 2)
b.Sigma <- rbind(c(1.6, .3, .3),
                 c(.3, 0.3, .4),
                 c(.3, .4, 2))
b.Sigmae2 <- rbind(c(.2, 0, 0),
                   c(0, .3, 0),
                   c(0, 0, .4))

# First, specify the ALpha1,Alpha2, theta, Sigma and sigmae2 parameters for each regime.
# Then we use the abind function to stack the parameters into arrays which's first
# dimension is the regime


Alpha <- abind::abind(a.Alpha, a.Alpha, along=-1, new.names=list(regime=c('a','a2'), x=NULL, y=NULL))
Theta <- abind::abind(a.Theta, a.Theta, along=-1, new.names=list(regime=c('a', 'a2'), xy=NULL))
Sigma <- abind::abind(a.Sigma, a.Sigma, along=-1, new.names=list(regime=c('a','a2'), x=NULL, y=NULL))
Sigmae <- abind::abind(a.Sigmae2, a.Sigmae2, along=-1, new.names=list(regime=c('a','a2'), x=NULL, y=NULL))


# regime 'a', traits 1, 2 and 3
model.a.123 <- list(X0 = a.X0,
                    Alpha=Alpha['a',,,drop=FALSE],
                    Theta=Theta['a',,drop=FALSE],
                    Sigma=Sigma['a',,,drop=FALSE],
                    Sigmae=Sigmae['a',,,drop=FALSE])
class(model.a.123) <- 'OU'

model.a.123.2SpOU <- list(X0 = a.X0,
                    Alpha1=Alpha['a2',,,drop=FALSE],
                    Alpha2=Alpha['a2',,,drop=FALSE],
                    Theta=Theta['a2',,drop=FALSE],
                    Sigma=Sigma['a2',,,drop=FALSE],
                    Sigmae=Sigmae['a2',,,drop=FALSE])
class(model.a.123.2SpOU) <- '2SpOU'

#####################################################################################################

context(ctx <- "R=1/k=3/N=400")

# number of tips
N <- 400

# tree with one regime
tree.a <- phytools::pbtree(n=N, scale=1)

# generate traits

traits.a.123 <- mvsim(tree.a, model.a.123, c(0,0,0), verbose=TRUE)

## Calculate likelihood

lik.OU = mvlik(traits.a.123$values+traits.a.123$errors, tree.a, model.a.123)

lik.2SpOU = mvlik(traits.a.123$values+traits.a.123$errors, tree.a, model.a.123.2SpOU)

cat('OU likelihood=',lik.OU,'\n')
cat('2SpOU likelihood=',lik.2SpOU,'\n')

test_that(paste(ctx, "Match multivariate likelihood of independent traits regime a"), {
  expect_true(abs(lik.OU - lik.2SpOU) < EPS)
})

#####################################################################################################

Alpha <- abind::abind(b.Alpha, b.Alpha, along=-1, new.names=list(regime=c('b','b2'), x=NULL, y=NULL))
Theta <- abind::abind(b.Theta, b.Theta, along=-1, new.names=list(regime=c('b', 'b2'), xy=NULL))
Sigma <- abind::abind(b.Sigma, b.Sigma, along=-1, new.names=list(regime=c('b','b2'), x=NULL, y=NULL))
Sigmae <- abind::abind(b.Sigmae2, b.Sigmae2, along=-1, new.names=list(regime=c('b','b2'), x=NULL, y=NULL))

# regime 'a', traits 1, 2 and 3
model.b.123 <- list(X0 = b.X0,
                    Alpha=Alpha['b',,,drop=FALSE],
                    Theta=Theta['b',,drop=FALSE],
                    Sigma=Sigma['b',,,drop=FALSE],
                    Sigmae=Sigmae['b',,,drop=FALSE])
class(model.b.123) <- 'OU'

model.b.123.2SpOU <- list(X0 = b.X0,
                          Alpha1=Alpha['b2',,,drop=FALSE],
                          Alpha2=Alpha['b2',,,drop=FALSE],
                          Theta=Theta['b2',,drop=FALSE],
                          Sigma=Sigma['b2',,,drop=FALSE],
                          Sigmae=Sigmae['b2',,,drop=FALSE])
class(model.b.123.2SpOU) <- '2SpOU'

context(ctx <- "R=1/k=3/N=400")

# number of tips
N <- 400

# tree with one regime
tree.b <- phytools::pbtree(n=N, scale=1)

# generate traits

traits.b.123 <- mvsim(tree.b, model.b.123, c(0,0,0), verbose=TRUE)

## Calculate likelihood

lik.OU = mvlik(traits.b.123$values+traits.b.123$errors, tree.b, model.b.123)

lik.2SpOU = mvlik(traits.b.123$values+traits.b.123$errors, tree.b, model.b.123.2SpOU)

cat('OU likelihood=',lik.OU,'\n')
cat('2SpOU likelihood=',lik.2SpOU,'\n')

test_that(paste(ctx, "Match multivariate likelihood of dependent traits regime b"), {
  expect_true(abs(lik.OU - lik.2SpOU) < EPS)
})
