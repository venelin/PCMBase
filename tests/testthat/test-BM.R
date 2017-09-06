library(phytools)
library(testthat)


set.seed(1)


#numerical zero:
EPS <- 10^-8

# number of regimes
R <- 2
# number of traits
k <- 3

# rate mtrix of transition from one regime to another
Q <- matrix(c(-1, 1, 1, -1), R, R)
colnames(Q) <- rownames(Q) <- letters[1:R]

## Specifying a bivariate BM process for each regime
# First, specify Sigma and sigmae2 parameters for each regime.
# Then we use the abind function to stack the parameters into arrays which's first
# dimension is the regime

# regimes



# in regime 'a' the three traits evolve according to three independent BM processes
a.X0 <- c(5, 2, 1)
a.Sigma <- rbind(c(1.6, 0, 0),
                 c(0, 2.4, 0),
                 c(0, 0, 2))
a.Sigmae2 <- rbind(c(0, 0, 0),
                   c(0, 0, 0),
                   c(0, 0, 0))

# in regime 'b' there is correlation between the traits. They evolve under the BM process.
b.X0 <- c(12, 4, 3)
b.Sigma <- rbind(c(1.6, .3, .3),
                 c(.3, 0.3, .4),
                 c(.3, .4, 2))
b.Sigmae2 <- rbind(c(.2, 0, 0),
                   c(0, .3, 0),
                   c(0, 0, .4))


a.OU.Alpha <- rbind(c(EPS, 0, 0),
                    c(0, EPS, 0),
                    c(0, 0, EPS))
a.OU.Theta <- c(0, 0, 0)

# in regime 'b.OU' there is correlation between the traits. They evolve under the OU process.
Alpha <- abind::abind(a.OU.Alpha, a.OU.Alpha, along=-1, new.names=list(regime=c('a','a.OU'), x=NULL, y=NULL))
Theta <- abind::abind(a.OU.Theta, a.OU.Theta, along=-1, new.names=list(regime=c('a', 'a.OU'), xy=NULL))
Sigma <- abind::abind(a.Sigma, a.Sigma, along=-1, new.names=list(regime=c('a','a.OU'), x=NULL, y=NULL))
Sigmae <- abind::abind(a.Sigmae2, a.Sigmae2, along=-1, new.names=list(regime=c('a','a.OU'), x=NULL, y=NULL))

#Sigma <- abind::abind(a.Sigma, b.Sigma, along=-1, new.names=list(regime=c('a','b'), x=NULL, y=NULL))
#Sigmae <- abind::abind(a.Sigmae2, b.Sigmae2, along=-1, new.names=list(regime=c('a','b'), x=NULL, y=NULL))

# regime 'a', traits 1, 2 and 3

# regime 'a', trait 1
model.a.1 <- list(X0 = a.X0[1],
                  Sigma=Sigma['a',1,1,drop=FALSE],
                  Sigmae=Sigmae['a',1,1,drop=FALSE])
class(model.a.1) <- 'BM'

model.a.1.OU <- list(X0 = a.X0[1],
                  Alpha=Alpha['a.OU',1,1,drop=FALSE],
                  Theta=Theta['a.OU',1,drop=FALSE],
                  Sigma=Sigma['a.OU',1,1,drop=FALSE],
                  Sigmae=Sigmae['a.OU',1,1,drop=FALSE])
class(model.a.1.OU) <- 'OU'


model.a.123 <- list(X0 = a.X0,
                    Sigma=Sigma['a',,,drop=FALSE],
                    Sigmae=Sigmae['a',,,drop=FALSE])
class(model.a.123) <- 'BM'

#print (dimnames(model.b.123.OU$Alpha)[[1]])

################ 1st Validation ######################################################



########
N <- 20

# tree with one regime
tree.a <- phytools::pbtree(n=N, scale=1)

traits.a.1 <- mvsim(tree.a, model.a.1, 0, verbose=TRUE)

mvlik(traits.a.1$values + traits.a.1$errors, tree.a,model.a.1)

mvlik(traits.a.1$values + traits.a.1$errors, tree.a,model.a.1.OU)

POUMM::likPOUMMGivenTreeVTips(
  z = traits.a.1$values[,1] + traits.a.1$errors[,1],
  tree = tree.a,
  0,
  0,
  sigma = sqrt(model.a.1$Sigma[1,1,1]),
  sigmae = sqrt(model.a.1$Sigmae[1,1,1]),
  g0 = a.X0[1])

#######

context(ctx <- "R=1/k=3/N=2")

# number of tips
N <- 20

# tree with one regime
tree.a <- phytools::pbtree(n=N, scale=1)

# generate traits

traits.a.123 <- mvsim(tree.a, model.a.123, c(0,0,0), verbose=TRUE)

BMlik = mvlik(X = traits.a.123$values + traits.a.123$errors,
      tree = tree.a,
      model = model.a.123)

POUMMlik = (POUMM::likPOUMMGivenTreeVTips(
  traits.a.123$values[,1]+traits.a.123$errors[,1],
  tree.a,
  0,
  0,
  sqrt(model.a.123$Sigma[1,1,1]),
  sqrt(model.a.123$Sigmae[1,1,1]),
  a.X0[1]) +

  POUMM::likPOUMMGivenTreeVTips(
    traits.a.123$values[,2]+traits.a.123$errors[,2],
    tree.a,
    0,
    0,
    sqrt(model.a.123$Sigma[1,2,2]),
    sqrt(model.a.123$Sigmae[1,2,2]),
    a.X0[2]) +

  POUMM::likPOUMMGivenTreeVTips(traits.a.123$values[,3]+traits.a.123$errors[,3],
                                tree.a,
                                0,
                                0,
                                sqrt(model.a.123$Sigma[1,3,3]),
                                sqrt(model.a.123$Sigmae[1,3,3]),
                                a.X0[3]))

cat('BM likelihood=',BMlik,'\n')
cat('POUMM likelihood=',POUMMlik,'\n')

test_that(paste(ctx, "Match multivariate likelihood of independent traits regime a"), {
  expect_true(abs(BMlik - POUMMlik) < EPS)
})


############################# Change Number of tips

context(ctx <- "R=1/k=3/N=400")

# number of tips
N <- 100

tree.a <- phytools::pbtree(n=N, scale=1)

# generate traits

traits.a.123 <- mvsim(tree.a, model.a.123, c(0,0,0), verbose=TRUE)

# tree with one regime

BMlik = mvlik(X = traits.a.123$values + traits.a.123$errors,
              tree = tree.a,
              model = model.a.123)

POUMMlik = (POUMM::likPOUMMGivenTreeVTips(
  traits.a.123$values[,1]+traits.a.123$errors[,1],
  tree.a,
  0,
  0,
  sqrt(model.a.123$Sigma[1,1,1]),
  sqrt(model.a.123$Sigmae[1,1,1]),
  a.X0[1]) +

    POUMM::likPOUMMGivenTreeVTips(
      traits.a.123$values[,2]+traits.a.123$errors[,2],
      tree.a,
      0,
      0,
      sqrt(model.a.123$Sigma[1,2,2]),
      sqrt(model.a.123$Sigmae[1,2,2]),
      a.X0[2]) +

    POUMM::likPOUMMGivenTreeVTips(traits.a.123$values[,3]+traits.a.123$errors[,3],
                                  tree.a,
                                  0,
                                  0,
                                  sqrt(model.a.123$Sigma[1,3,3]),
                                  sqrt(model.a.123$Sigmae[1,3,3]),
                                  a.X0[3]))

cat('BM likelihood=',BMlik,'\n')
cat('POUMM likelihood=',POUMMlik,'\n')

test_that(paste(ctx, "Match multivariate likelihood of independent traits regime a"), {
  expect_true(abs(BMlik - POUMMlik) < EPS)
})


######################################################################################

################ 2nd Validation ######################################################

EPS <- 10^-8

b.X0 <- c(12, 4, 3)
b.Sigma <- rbind(c(1.6, .3, .3),
                 c(.3, 0.3, .4),
                 c(.3, .4, 2))
b.Sigmae2 <- rbind(c(0, 0, 0),
                   c(0, 0, 0),
                   c(0, 0, 0))

b.OU.Alpha <- rbind(c(EPS, 0, 0),
                    c(0, EPS, 0),
                    c(0, 0, EPS))
b.OU.Theta <- c(0, 0, 0)


Alpha <- abind::abind(b.OU.Alpha, b.OU.Alpha, along=-1, new.names=list(regime=c('b','b.OU'), x=NULL, y=NULL))
Theta <- abind::abind(b.OU.Theta, b.OU.Theta, along=-1, new.names=list(regime=c('b', 'b.OU'), xy=NULL))
Sigma <- abind::abind(b.Sigma, b.Sigma, along=-1, new.names=list(regime=c('b','b.OU'), x=NULL, y=NULL))
Sigmae <- abind::abind(b.Sigmae2, b.Sigmae2, along=-1, new.names=list(regime=c('b','b.OU'), x=NULL, y=NULL))

model.b.123 <- list(X0 = b.X0,
                    Sigma=Sigma['b',,,drop=FALSE],
                    Sigmae=Sigmae['b',,,drop=FALSE])
class(model.b.123) <- 'BM'

model.b.123.OU <- list(X0 = b.X0,
                       Alpha=Alpha['b.OU',,,drop=FALSE],
                       Theta=Theta['b.OU',,drop=FALSE],
                       Sigma=Sigma['b.OU',,,drop=FALSE],
                       Sigmae=Sigmae['b.OU',,,drop=FALSE])
class(model.b.123.OU) <- 'OU'

context(ctx <- "R=1/k=3/N=2")

N=10

tree.b <- phytools::pbtree(n=N, scale=1)

traits.b.123 <- mvsim(tree.b, model.b.123, c(0,0,0), verbose=TRUE)

## Calculate likelihood

lik.BM = mvlik(X = traits.b.123$values + traits.b.123$errors, tree = tree.b, model = model.b.123)

lik.OU = mvlik(X = traits.b.123$values + traits.b.123$errors, tree = tree.b, model = model.b.123.OU)


POUMMlik = (POUMM::likPOUMMGivenTreeVTips(
  traits.b.123$values[,1]+traits.b.123$errors[,1],
  tree.b,
  0,
  0,
  sqrt(model.b.123$Sigma[1,1,1]),
  sqrt(model.b.123$Sigmae[1,1,1]),
  b.X0[1]) +

    POUMM::likPOUMMGivenTreeVTips(
      traits.b.123$values[,2]+traits.b.123$errors[,2],
      tree.b,
      0,
      0,
      sqrt(model.b.123$Sigma[1,2,2]),
      sqrt(model.b.123$Sigmae[1,2,2]),
      b.X0[2]) +

    POUMM::likPOUMMGivenTreeVTips(traits.b.123$values[,3]+traits.b.123$errors[,3],
                                  tree.b,
                                  0,
                                  0,
                                  sqrt(model.b.123$Sigma[1,3,3]),
                                  sqrt(model.b.123$Sigmae[1,3,3]),
                                  b.X0[3]))


cat('BM likelihood=',lik.BM,'\n')
cat('OU likelihood=',lik.OU,'\n')
cat('POUMM likelihood=',POUMMlik,'\n')


test_that(paste(ctx, "Match multivariate likelihood of independent traits regime b"),{
  expect_true( abs(lik.BM - lik.OU ) < EPS)
  })

#############################################################################

tree.ab <- phytools::sim.history(tree.a, Q, anc='a')

# convert the simmap tree to a normal phylo object with singleton nodes at the
# within-branch regime changes. The regimes are encoded as names of the edge.length
# vector
tree.ab.singles <- map.to.singleton(tree.ab)


traits.ab.123 <- mvsim(tree.ab, model.ab.123, c(0,0,0), verbose=TRUE)

