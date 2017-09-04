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

Alpha <- abind::abind(a.Alpha, b.Alpha, along=-1, new.names=list(regime=c('a','b'), x=NULL, y=NULL))
Theta <- abind::abind(a.Theta, b.Theta, along=-1, new.names=list(regime=c('a', 'b'), xy=NULL))
Sigma <- abind::abind(a.Sigma, b.Sigma, along=-1, new.names=list(regime=c('a','b'), x=NULL, y=NULL))
Sigmae <- abind::abind(a.Sigmae2, b.Sigmae2, along=-1, new.names=list(regime=c('a','b'), x=NULL, y=NULL))


## Simulations of trait data

# regime 'a', trait 1
model.a.1 <- list(X0 = a.X0[1],
                  Alpha=Alpha['a',1,1,drop=FALSE],
                  Theta=Theta['a',1,drop=FALSE],
                  Sigma=Sigma['a',1,1,drop=FALSE],
                  Sigmae=Sigmae['a',1,1,drop=FALSE])
class(model.a.1) <- 'OU'

# regime 'a', trait 2
model.a.2 <- list(X0 = a.X0[2],
                  Alpha=Alpha['a',2,2,drop=FALSE],
                  Theta=Theta['a',2,drop=FALSE],
                  Sigma=Sigma['a',2,2,drop=FALSE],
                  Sigmae=Sigmae['a',2,2,drop=FALSE])
class(model.a.2) <- 'OU'

# regime 'a', trait 3
model.a.3 <- list(X0 = a.X0[3],
                  Alpha=Alpha['a',3,3,drop=FALSE],
                  Theta=Theta['a',3,drop=FALSE],
                  Sigma=Sigma['a',3,3,drop=FALSE],
                  Sigmae=Sigmae['a',3,3,drop=FALSE])
class(model.a.3) <- 'OU'

# regime 'a', traits 1, 2 and 3
model.a.123 <- list(X0 = a.X0,
                    Alpha=Alpha['a',,,drop=FALSE],
                    Theta=Theta['a',,drop=FALSE],
                    Sigma=Sigma['a',,,drop=FALSE],
                    Sigmae=Sigmae['a',,,drop=FALSE])
class(model.a.123) <- 'OU'


# regime 'b', traits 1, 2 and 3
model.b.123 <- list(X0 = b.X0,
                    Alpha=Alpha['b',,,drop=FALSE],
                    Theta=Theta['b',,drop=FALSE],
                    Sigma=Sigma['b',,,drop=FALSE],
                    Sigmae=Sigmae['b',,,drop=FALSE])
class(model.b.123) <- 'OU'

# regimes 'a' and 'b', traits 1, 2 and 3
model.ab.123 <- list(X0 = a.X0,
                     Alpha=Alpha[,,,drop=FALSE],
                     Theta=Theta[,,drop=FALSE],
                     Sigma=Sigma[,,,drop=FALSE],
                     Sigmae=Sigmae[,,,drop=FALSE])
class(model.ab.123) <- 'OU'


context(ctx <- "R=1/k=1/N=2")

# number of tips
N <- 2

# tree with one regime
tree.a <- phytools::pbtree(n=N, scale=1)

# generate traits
traits.a.1 <- mvsim(tree.a, model.a.1, 0, verbose=TRUE)
traits.a.2 <- mvsim(tree.a, model.a.2, 0, verbose=TRUE)
traits.a.3 <- mvsim(tree.a, model.a.2, 0, verbose=TRUE)


# test likelihood
test_that(paste(ctx, "Match univariate likelihood from patherit regime a"), {
  expect_true(
    abs(mvlik(traits.a.1$values + traits.a.1$errors, tree.a,
            model.a.1) -
          POUMM::likPOUMMGivenTreeVTips(
            z = traits.a.1$values[,1] + traits.a.1$errors[,1],
            tree = tree.a,
            alpha = model.a.1$A[1,1,1],
            theta = model.a.1$Theta[1,1],
            sigma = sqrt(model.a.1$Sigma[1,1,1]),
            sigmae = sqrt(model.a.1$Sigmae[1,1,1]),
            g0 = a.X0[1])
        ) < EPS)
  expect_true(
    abs(mvlik(traits.a.2$values + traits.a.2$errors, tree.a,
        model.a.2) -
          POUMM::likPOUMMGivenTreeVTips(
            z = traits.a.2$values[,1] + traits.a.2$errors[,1],
            tree = tree.a,
            alpha = model.a.2$A[1,1,1],
            theta = model.a.2$Theta[1,1],
            sigma = sqrt(model.a.2$Sigma[1,1,1]),
            sigmae = sqrt(model.a.2$Sigmae[1,1,1]),
            g0 = a.X0[2])
        ) < EPS)
  expect_true(
    abs(mvlik(traits.a.3$values + traits.a.3$errors, tree.a,
        model.a.3) -
          POUMM::likPOUMMGivenTreeVTips(
            z = traits.a.3$values[,1] + traits.a.3$errors[,1],
            tree = tree.a,
            alpha = model.a.3$A[1,1,1],
            theta = model.a.3$Theta[1,1],
            sigma = sqrt(model.a.3$Sigma[1,1,1]),
            sigmae = sqrt(model.a.3$Sigmae[1,1,1]),
            g0 = a.X0[3])) < EPS)
})

context(ctx <- "R=1/k=1/N=400")

# number of tips
N <- 400

# tree with one regime
tree.a <- pbtree(n=N, scale=1)

# generate traits
traits.a.1 <- mvsim(tree.a, model.a.1, 0, verbose=TRUE)
traits.a.2 <- mvsim(tree.a, model.a.2, 0, verbose=TRUE)
traits.a.3 <- mvsim(tree.a, model.a.2, 0, verbose=TRUE)

test_that(paste(ctx, "Match univariate likelihood from patherit regime a"), {
  expect_true(
    abs(mvlik(traits.a.1$values + traits.a.1$errors, tree.a,
              model.a.1) -
          POUMM::likPOUMMGivenTreeVTips(
            z = traits.a.1$values[,1] + traits.a.1$errors[,1],
            tree = tree.a,
            alpha = model.a.1$A[1,1,1],
            theta = model.a.1$Theta[1,1],
            sigma = sqrt(model.a.1$Sigma[1,1,1]),
            sigmae = sqrt(model.a.1$Sigmae[1,1,1]),
            g0 = a.X0[1])
    ) < EPS)
  expect_true(
    abs(mvlik(traits.a.2$values + traits.a.2$errors, tree.a,
              model.a.2) -
          POUMM::likPOUMMGivenTreeVTips(
            z = traits.a.2$values[,1] + traits.a.2$errors[,1],
            tree = tree.a,
            alpha = model.a.2$A[1,1,1],
            theta = model.a.2$Theta[1,1],
            sigma = sqrt(model.a.2$Sigma[1,1,1]),
            sigmae = sqrt(model.a.2$Sigmae[1,1,1]),
            g0 = a.X0[2])
    ) < EPS)
  expect_true(
    abs(mvlik(traits.a.3$values + traits.a.3$errors, tree.a,
              model.a.3) -
          POUMM::likPOUMMGivenTreeVTips(
            z = traits.a.3$values[,1] + traits.a.3$errors[,1],
            tree = tree.a,
            alpha = model.a.3$A[1,1,1],
            theta = model.a.3$Theta[1,1],
            sigma = sqrt(model.a.3$Sigma[1,1,1]),
            sigmae = sqrt(model.a.3$Sigmae[1,1,1]),
            g0 = a.X0[3])) < EPS)
})

context(ctx <- "R=1/k=3/N=2")

traits.a.123 <- mvsim(tree.a, model.a.123, c(0,0,0), verbose=TRUE)

test_that(paste(ctx, "Match multivariate likelihood of independent traits regime a"), {
  expect_true(
    abs(mvlik(traits.a.123$values+traits.a.123$errors, tree.a, model.a.123) -
          (POUMM::likPOUMMGivenTreeVTips(
            traits.a.123$values[,1]+traits.a.123$errors[,1],
            tree.a,
            model.a.123$A[1,1,1],
            model.a.123$Theta[1,1],
            sqrt(model.a.123$Sigma[1,1,1]),
            sqrt(model.a.123$Sigmae[1,1,1]),
            a.X0[1]) +

             POUMM::likPOUMMGivenTreeVTips(
               traits.a.123$values[,2]+traits.a.123$errors[,2],
               tree.a,
               model.a.123$A[1,2,2],
               model.a.123$Theta[1,2],
               sqrt(model.a.123$Sigma[1,2,2]),
               sqrt(model.a.123$Sigmae[1,2,2]),
               a.X0[2]) +

             POUMM::likPOUMMGivenTreeVTips(traits.a.123$values[,3]+traits.a.123$errors[,3],
                              tree.a,
                              model.a.123$A[1,3,3],
                              model.a.123$Theta[1,3],
                              sqrt(model.a.123$Sigma[1,3,3]),
                              sqrt(model.a.123$Sigmae[1,3,3]),
                              a.X0[3]))
        ) < EPS)
})



context(ctx <- "R=1/k=3/N=400")

traits.a.123 <- mvsim(tree.a, model.a.123, c(0,0,0), verbose=TRUE)


## Calculate likelihood
test_that(paste(ctx, "Match multivariate likelihood of independent traits regime a"), {
  expect_true(
    abs(mvlik(traits.a.123$values+traits.a.123$errors, tree.a, model.a.123) -
          (POUMM::likPOUMMGivenTreeVTips(
            traits.a.123$values[,1]+traits.a.123$errors[,1],
            tree.a,
            model.a.123$A[1,1,1],
            model.a.123$Theta[1,1],
            sqrt(model.a.123$Sigma[1,1,1]),
            sqrt(model.a.123$Sigmae[1,1,1]),
            a.X0[1]) +
             POUMM::likPOUMMGivenTreeVTips(
               traits.a.123$values[,2] + traits.a.123$errors[,2],
               tree.a,
               model.a.123$A[1,2,2],
               model.a.123$Theta[1,2],
               sqrt(model.a.123$Sigma[1,2,2]),
               sqrt(model.a.123$Sigmae[1,2,2]),
               a.X0[2]) +
             POUMM::likPOUMMGivenTreeVTips(
               traits.a.123$values[,3]+traits.a.123$errors[,3],
               tree.a,
               model.a.123$A[1,3,3],
               model.a.123$Theta[1,3],
               sqrt(model.a.123$Sigma[1,3,3]),
               sqrt(model.a.123$Sigmae[1,3,3]),
               a.X0[3]))
        ) < EPS)
})


plot(traits.a.1$values)


traits.b.123 <- mvsim(tree.a, model.b.123, c(0,0,0), verbose=TRUE)

tree.ab <- phytools::sim.history(tree.a, Q, anc='a')

# convert the simmap tree to a normal phylo object with singleton nodes at the
# within-branch regime changes. The regimes are encoded as names of the edge.length
# vector
tree.ab.singles <- map.to.singleton(tree.ab)


traits.ab.123 <- mvsim(tree.ab, model.ab.123, c(0,0,0), verbose=TRUE)


