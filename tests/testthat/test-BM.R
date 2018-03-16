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

model.a.123 <- PCM("BM", k = 3, regimes = "a",
                   params = list(X0 = a.X0, Sigma = Sigma[,, "a", drop = FALSE],
                                 Sigmae = Sigmae[,, "a", drop = FALSE]))


################ 1st Validation ######################################################

context(ctx <- "R=1/k=3/N=2")

# number of tips
N <- 2

# tree with one regime
tree.a <- rtree(N) # phytools::pbtree(n=N, scale=1)
PCMSetDefaultRegime(tree.a, model.a.123)

# generate traits

traits.a.123 <- PCMSim(tree.a, model.a.123, c(0,0,0), verbose=TRUE)

BMlik = PCMLik(X = traits.a.123[,1:N],
      tree = tree.a,
      model = model.a.123)

POUMMlik = (POUMM::likPOUMMGivenTreeVTips(
  traits.a.123[1,],
  tree.a,
  0,
  0,
  sqrt(model.a.123$Sigma[1,1,1]),
  sqrt(model.a.123$Sigmae[1,1,1]),
  a.X0[1]) +

  POUMM::likPOUMMGivenTreeVTips(
    traits.a.123[2,],
    tree.a,
    0,
    0,
    sqrt(model.a.123$Sigma[2,2,1]),
    sqrt(model.a.123$Sigmae[2,2,1]),
    a.X0[2]) +

  POUMM::likPOUMMGivenTreeVTips(traits.a.123[3,],
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
PCMSetDefaultRegime(tree.a, model.a.123)
# generate traits

traits.a.123 <- PCMSim(tree.a, model.a.123, c(0,0,0), verbose=TRUE)

# tree with one regime

BMlik = PCMLik(X = traits.a.123,
              tree = tree.a,
              model = model.a.123)

POUMMlik = (POUMM::likPOUMMGivenTreeVTips(
  traits.a.123[1,],
  tree.a,
  0,
  0,
  sqrt(model.a.123$Sigma[1,1,1]),
  sqrt(model.a.123$Sigmae[1,1,1]),
  a.X0[1]) +

    POUMM::likPOUMMGivenTreeVTips(
      traits.a.123[2,],
      tree.a,
      0,
      0,
      sqrt(model.a.123$Sigma[2,2,1]),
      sqrt(model.a.123$Sigmae[2,2,1]),
      a.X0[2]) +

    POUMM::likPOUMMGivenTreeVTips(traits.a.123[3,],
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

model.b.123 <- PCM("BM", k = 3, regimes = "b",
                   params = list(X0 = b.X0,
                                 Sigma=Sigma[,,'b',drop=FALSE],
                                 Sigmae=Sigmae[,,'b',drop=FALSE]))


model.b.123.OU <- PCM("OU", k = 3, regimes = "b", params =
                        list(X0 = b.X0,
                             H=H[,,'b.OU',drop=FALSE],
                             Theta=Theta[,'b.OU',drop=FALSE],
                             Sigma=Sigma[,,'b.OU',drop=FALSE],
                             Sigmae=Sigmae[,,'b.OU',drop=FALSE]))


context(ctx <- "R=1/k=3/N=400")

N=400

tree.b <- rtree(N) # phytools::pbtree(n=N, scale=1)

PCMSetDefaultRegime(tree.b, model.b.123)

traits.b.123 <- PCMSim(tree.b, model.b.123, c(0,0,0), verbose=TRUE)

## Calculate likelihood
lik.BM = PCMLik(X = traits.b.123, tree = tree.b, model = model.b.123)
lik.OU = PCMLik(X = traits.b.123, tree = tree.b, model = model.b.123.OU)

POUMMlik = (POUMM::likPOUMMGivenTreeVTips(
  traits.b.123[1,],
  tree.b,
  0,
  0,
  sqrt(model.b.123$Sigma[1,1,1]),
  sqrt(model.b.123$Sigmae[1,1,1]),
  b.X0[1]) +

    POUMM::likPOUMMGivenTreeVTips(
      traits.b.123[2,],
      tree.b,
      0,
      0,
      sqrt(model.b.123$Sigma[2,2,1]),
      sqrt(model.b.123$Sigmae[2,2,1]),
      b.X0[2]) +

    POUMM::likPOUMMGivenTreeVTips(traits.b.123[3,],
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

### WHITE MODEL TEST ###
model.a.123.white <- PCMSetParams(model.a.123, params = list(Sigma = abind(matrix(0, 3, 3), along = 3),
                                        Sigmae = abind(diag(1:3, 3, 3), along = 3)), inplace = FALSE)

lik.white <- PCMLik(traits.a.123[,1:N], tree.a, model.a.123.white)
cat("White likelihood=", lik.white, "\n")

lik.white.true <- sum(
  apply(traits.a.123[, 1:N], 2,
        function(xi) mvtnorm::dmvnorm(x = xi, mean = c(5, 2, 1),
                                      sigma = model.a.123.white$Sigmae[,,1], log = TRUE )))

test_that("Match white likelihoods", expect_true(abs(lik.white-lik.white.true) < EPS))


if(require(PCMBaseCpp)) {
  cat("Testing PCMBaseCpp on BM:\n")

  lik.white <- PCMLik(tree = tree.a, model = model.a.123.white,
                      metaI = PCMInfoCpp(X = traits.a.123[, 1:N], tree = tree.a, model = model.a.123.white))

  test_that("Match white likelihoods", expect_true(abs(lik.white-lik.white.true) < EPS))

  test_that("a.123",
            expect_equal(PCMLik(traits.a.123, tree.a, model.a.123),
                         PCMLik(traits.a.123, tree.a, model.a.123,
                               metaI = PCMInfoCpp(X = traits.a.123[, 1:length(tree.a$tip.label)],
                                                     tree = tree.a,
                                                     model.a.123))))

  values <- traits.a.123[, 1:length(tree.a$tip.label)]
  values[sample(x=1:length(values), 50)] <- NA

  metaIR <- PCMInfo(X = values, tree = tree.a, model = model.a.123)
  metaICpp <- PCMInfoCpp(X = values, tree = tree.a, model = model.a.123)

  test_that("a.123 with missing values",
            expect_equal(PCMLik(values, tree.a, model.a.123),
                         PCMLik(tree = tree.a, model = model.a.123, metaI = metaICpp)))




  if(require(microbenchmark)) {
    cat("microbenchmark test")

    options(PCMBase.PCMLmr.mode=11)
    print(microbenchmark(
      PCMLik(values, tree.a, model.a.123, metaI = metaIR),
      PCMLik(values, tree.a, model.a.123, metaI = metaICpp), times = 10
    ))

    options(PCMBase.PCMLmr.mode=21)
    print(microbenchmark(
      PCMLik(values, tree.a, model.a.123, metaI = metaICpp)
    ))
  }

}

