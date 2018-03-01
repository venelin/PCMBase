library(ape)
library(testthat)
library(PCMBase)
library(abind)

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

# First, specify Sigma and sigmae2 parameters for each regime.
# Then we use the abind function to stack the parameters into arrays which's first
# dimension is the regime

Sigma <- abind(a.Sigma, b.Sigma, along=3, new.names=list(x=NULL, y=NULL, regime=c('a','b')))
Sigmae <- abind(a.Sigmae2, b.Sigmae2, along=3, new.names=list(x=NULL, y=NULL, regime=c('a','b')))

# regime 'a', traits 1, 2 and 3

model.a.123 <- list(X0 = a.X0,
                    Sigma=Sigma[,,'a',drop=FALSE],
                    Sigmae=Sigmae[,,'a',drop=FALSE])
class(model.a.123) <- 'BM'

################ 1st Validation ######################################################

context(ctx <- "R=1/k=3/N=2")

# number of tips
N <- 2

# tree with one regime
tree.a <- rtree(N) # phytools::pbtree(n=N, scale=1)

# generate traits

traits.a.123 <- PCMSim(tree.a, model.a.123, c(0,0,0), verbose=TRUE)

BMlik = PCMLik(X = traits.a.123$values + traits.a.123$errors,
      tree = tree.a,
      model = model.a.123)

POUMMlik = (POUMM::likPOUMMGivenTreeVTips(
  traits.a.123$values[1,]+traits.a.123$errors[1,],
  tree.a,
  0,
  0,
  sqrt(model.a.123$Sigma[1,1,1]),
  sqrt(model.a.123$Sigmae[1,1,1]),
  a.X0[1]) +

  POUMM::likPOUMMGivenTreeVTips(
    traits.a.123$values[2,]+traits.a.123$errors[2,],
    tree.a,
    0,
    0,
    sqrt(model.a.123$Sigma[2,2,1]),
    sqrt(model.a.123$Sigmae[2,2,1]),
    a.X0[2]) +

  POUMM::likPOUMMGivenTreeVTips(traits.a.123$values[3,]+traits.a.123$errors[3,],
                                tree.a,
                                0,
                                0,
                                sqrt(model.a.123$Sigma[3,3,1]),
                                sqrt(model.a.123$Sigmae[3,3,1]),
                                a.X0[3]))

cat('BM likelihood=',BMlik,'\n')
cat('POUMM likelihood=',POUMMlik,'\n')

test_that(paste(ctx, "Match multivariate likelihood of independent traits regime a"), {
  expect_true(abs(BMlik - POUMMlik) < EPS)
})


############################# Change Number of tips

context(ctx <- "R=1/k=3/N=400")

# number of tips
N <- 400

tree.a <- rtree(N) # phytools::pbtree(n=N, scale=1)

# generate traits

traits.a.123 <- PCMSim(tree.a, model.a.123, c(0,0,0), verbose=TRUE)

# tree with one regime

BMlik = PCMLik(X = traits.a.123$values + traits.a.123$errors,
              tree = tree.a,
              model = model.a.123)

POUMMlik = (POUMM::likPOUMMGivenTreeVTips(
  traits.a.123$values[1,]+traits.a.123$errors[1,],
  tree.a,
  0,
  0,
  sqrt(model.a.123$Sigma[1,1,1]),
  sqrt(model.a.123$Sigmae[1,1,1]),
  a.X0[1]) +

    POUMM::likPOUMMGivenTreeVTips(
      traits.a.123$values[2,]+traits.a.123$errors[2,],
      tree.a,
      0,
      0,
      sqrt(model.a.123$Sigma[2,2,1]),
      sqrt(model.a.123$Sigmae[2,2,1]),
      a.X0[2]) +

    POUMM::likPOUMMGivenTreeVTips(traits.a.123$values[3,]+traits.a.123$errors[3,],
                                  tree.a,
                                  0,
                                  0,
                                  sqrt(model.a.123$Sigma[3,3,1]),
                                  sqrt(model.a.123$Sigmae[3,3,1]),
                                  a.X0[3]))

cat('BM likelihood=',BMlik,'\n')
cat('POUMM likelihood=',POUMMlik,'\n')

test_that(paste(ctx, "Match multivariate likelihood of independent traits regime a"), {
  expect_true(abs(BMlik - POUMMlik) < EPS)
})


######################################################################################

################ 2nd Validation ######################################################

b.OU.H <- rbind(c(EPS, 0, 0),
                    c(0, EPS, 0),
                    c(0, 0, EPS))
b.OU.Theta <- c(0, 0, 0)


H <- abind(b.OU.H, b.OU.H, along=3, new.names=list(x=NULL, y=NULL, regime=c('b','b.OU')))
Theta <- abind(b.OU.Theta, b.OU.Theta, along=2, new.names=list(xy=NULL, regime=c('b', 'b.OU')))
Sigma <- abind(b.Sigma, b.Sigma, along=3, new.names=list(x=NULL, y=NULL, regime=c('b','b.OU')))
Sigmae <- abind(b.Sigmae2, b.Sigmae2, along=3, new.names=list(x=NULL, y=NULL, regime=c('b','b.OU')))

model.b.123 <- list(X0 = b.X0,
                    Sigma=Sigma[,,'b',drop=FALSE],
                    Sigmae=Sigmae[,,'b',drop=FALSE])
class(model.b.123) <- 'BM'

model.b.123.OU <- list(X0 = b.X0,
                       H=H[,,'b.OU',drop=FALSE],
                       Theta=Theta[,'b.OU',drop=FALSE],
                       Sigma=Sigma[,,'b.OU',drop=FALSE],
                       Sigmae=Sigmae[,,'b.OU',drop=FALSE])
class(model.b.123.OU) <- 'OU'

context(ctx <- "R=1/k=3/N=400")

N=400

tree.b <- rtree(N) # phytools::pbtree(n=N, scale=1)

traits.b.123 <- PCMSim(tree.b, model.b.123, c(0,0,0), verbose=TRUE)

## Calculate likelihood
lik.BM = PCMLik(X = traits.b.123$values + traits.b.123$errors, tree = tree.b, model = model.b.123)
lik.OU = PCMLik(X = traits.b.123$values + traits.b.123$errors, tree = tree.b, model = model.b.123.OU)

POUMMlik = (POUMM::likPOUMMGivenTreeVTips(
  traits.b.123$values[1,]+traits.b.123$errors[1,],
  tree.b,
  0,
  0,
  sqrt(model.b.123$Sigma[1,1,1]),
  sqrt(model.b.123$Sigmae[1,1,1]),
  b.X0[1]) +

    POUMM::likPOUMMGivenTreeVTips(
      traits.b.123$values[2,]+traits.b.123$errors[2,],
      tree.b,
      0,
      0,
      sqrt(model.b.123$Sigma[2,2,1]),
      sqrt(model.b.123$Sigmae[2,2,1]),
      b.X0[2]) +

    POUMM::likPOUMMGivenTreeVTips(traits.b.123$values[3,]+traits.b.123$errors[3,],
                                  tree.b,
                                  0,
                                  0,
                                  sqrt(model.b.123$Sigma[3,3,1]),
                                  sqrt(model.b.123$Sigmae[3,3,1]),
                                  b.X0[3]))

cat('BM likelihood=',lik.BM,'\n')
cat('OU likelihood=',lik.OU,'\n')
cat('POUMM likelihood=',POUMMlik,'\n')


test_that(paste(ctx, "Match multivariate likelihood of independent traits regime b"),{
  expect_true( abs(lik.BM - lik.OU ) < EPS*1000)
  })

if(require(PCMBaseCpp)) {
  cat("Testing PCMBaseCpp on BM:\n")

  test_that("a.123",
            expect_equal(PCMLik(traits.a.123$values+traits.a.123$errors, tree.a, model.a.123),
                         PCMLik(traits.a.123$values+traits.a.123$errors, tree.a, model.a.123,
                               pruneI = PCMCppPruningObject(X = traits.a.123$values[, 1:length(tree.a$tip.label)],
                                                     tree = tree.a,
                                                     model.a.123))))

  values <- traits.a.123$values[, 1:length(tree.a$tip.label)] + traits.a.123$errors[, 1:length(tree.a$tip.label)]
  values[sample(x=1:length(values), 50)] <- NA

  pruneInfoR <- PCMPruningOrder(tree.a)
  pruneInfoCpp <- PCMCppPruningObject(X = values,
                               tree = tree.a,
                               model.a.123)

  test_that("a.123 with missing values",
            expect_equal(PCMLik(values, tree.a, model.a.123),
                         PCMLik(tree = tree.a, model = model.a.123,
                               pruneI = pruneInfoCpp)))


  # if(require(microbenchmark)) {
  #   cat("microbenchmark test")
  #
  #   options(PCMBase.PCMLmr.mode=11)
  #   print(microbenchmark(
  #     PCMLik(values, tree.a, model.a.123, pruneI = pruneInfoR),
  #     PCMLik(values, tree.a, model.a.123, pruneI = pruneInfoCpp), times = 10
  #   ))
  #
  #   options(PCMBase.PCMLmr.mode=21)
  #   print(microbenchmark(
  #     PCMLik(values, tree.a, model.a.123, pruneI = pruneInfoCpp)
  #   ))
  # }

}


