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

# regimes

# in regime 'a' the three traits evolve according to three independent OU processes
a.X0 <- c(5, 2, 1)
a.H <- rbind(
  c(0, 0, 0),
  c(0, 2, 0),
  c(0, 0, 3))
a.Theta <- c(10, 6, 2)
a.Sigma_x <- rbind(
  c(1.6, 0.0, 0.0),
  c(0.0, 2.4, 0.0),
  c(0.0, 0.0, 2.0))
a.Sigmae_x <- rbind(
  c(0.0, 0.0, 0.0),
  c(0.0, 0.0, 0.0),
  c(0.0, 0.0, 0.0))

# in regime 'b' there is correlation between the traits
b.X0 <- c(12, 4, 3)
b.H <- rbind(
  c(2.0, 0.1, 0.2),
  c(0.1, 0.6, 0.2),
  c(0.2, 0.2, 0.3))
b.Theta <- c(10, 6, 2)
b.Sigma_x <- rbind(
  c(1.6, 0.3, 0.3),
  c(0.0, 0.3, 0.4),
  c(0.0, 0.0, 2.0))
b.Sigmae_x <- rbind(
  c(0.2, 0.0, 0.0),
  c(0.0, 0.3, 0.0),
  c(0.0, 0.0, 0.4))

H <- abind(a.H, b.H, along=3, new.names=list(x=NULL, y=NULL, regime=c('a','b')))
Theta <- abind(a.Theta, b.Theta, along=2, new.names=list(xy=NULL, regime=c('a','b')))
Sigma_x <- abind(a.Sigma_x, b.Sigma_x, along=3, new.names=list(x=NULL, y=NULL, regime=c('a','b')))
Sigmae_x <- abind(a.Sigmae_x, b.Sigmae_x, along=3, new.names=list(x=NULL, y=NULL, regime=c('a','b')))


## Simulations of trait data

# regime 'a', trait 1
model.a.1 <- PCM("OU", k = 1, regimes = "a", params = list(X0 = a.X0[1],
                                             H=H[1,1,'a',drop=FALSE],
                                             Theta=Theta[1,'a',drop=FALSE],
                                             Sigma_x=Sigma_x[1,1,'a',drop=FALSE],
                                             Sigmae_x=Sigmae_x[1,1,'a',drop=FALSE]))


# regime 'a', trait 2
model.a.2 <- PCM("OU", k = 1, regimes = "a", params = list(X0 = a.X0[2],
                                             H=H[2,2,'a',drop=FALSE],
                                             Theta=Theta[2,'a',drop=FALSE],
                                             Sigma_x=Sigma_x[2,2,'a',drop=FALSE],
                                             Sigmae_x=Sigmae_x[2,2,'a',drop=FALSE]))

# regime 'a', trait 3
model.a.3 <- PCM("OU", k = 1, regimes = "a", params = list(X0 = a.X0[3],
                                             H=H[3,3,'a',drop=FALSE],
                                             Theta=Theta[3,'a',drop=FALSE],
                                             Sigma_x=Sigma_x[3,3,'a',drop=FALSE],
                                             Sigmae_x=Sigmae_x[3,3,'a',drop=FALSE]))

# regime 'a', traits 1, 2 and 3
model.a.123 <- PCM("OU", k = 3, regimes = "a", params = list(X0 = a.X0,
                                               H=H[,,'a',drop=FALSE],
                                               Theta=Theta[,'a',drop=FALSE],
                                               Sigma_x=Sigma_x[,,'a',drop=FALSE],
                                               Sigmae_x=Sigmae_x[,,'a',drop=FALSE]))


# regime 'b', traits 1, 2 and 3
model.b.123 <- PCM("OU", k = 3, regimes = "b", params = list(X0 = b.X0,
                                               H=H[,,'b',drop=FALSE],
                                               Theta=Theta[,'b',drop=FALSE],
                                               Sigma_x=Sigma_x[,,'b',drop=FALSE],
                                               Sigmae_x=Sigmae_x[,,'b',drop=FALSE]))

# regimes 'a' and 'b', traits 1, 2 and 3
model.ab.123 <- PCM("OU", k = 3, regimes = c("a", "b"), params = list(X0 = a.X0,
                                                        H=H[,,,drop=FALSE],
                                                        Theta=Theta[,,drop=FALSE],
                                                        Sigma_x=Sigma_x[,,,drop=FALSE],
                                                        Sigmae_x=Sigmae_x[,,,drop=FALSE]))

cat("PCMNumParams(model.ab.123)=", PCMNumParams(model.ab.123), "\n")
cat("length(PCMGetVecParamsFull(model.ab.123)=", length(PCMGetVecParamsFull(model.ab.123)), "\n")

test_that("Check correctness of PCMGetVecParamsFull",
          expect_equal(length(PCMGetVecParamsFull(model.ab.123)), 2*(3 + 3*3 + 3 + 3*3 + 3*3)))

test_that("Check correctness of PCMNumParams",
          expect_equal(PCMNumParams(model.ab.123), 3 + 2*(3*3 + 3 + 6 + 6)))


context(ctx <- "OU: R=1/k=1/N=2")

# number of tips
N <- 2

# tree with one regime
tree.a <- rtree(N) # phytools::pbtree(n=N, scale=1)
PCMTreeSetDefaultRegime(tree.a, model.a.123)
#tree.a$edge.regime <- rep("a", length(tree.a$edge.length))

# generate traits
traits.a.1 <- PCMSim(tree.a, model.a.1, 0, verbose=TRUE)
traits.a.2 <- PCMSim(tree.a, model.a.2, 0, verbose=TRUE)
traits.a.3 <- PCMSim(tree.a, model.a.2, 0, verbose=TRUE)


# test likelihood
test_that(paste(ctx, "Match univariate likelihood from patherit regime a"), {
  expect_true(
    abs(PCMLik(traits.a.1, tree.a,
            model.a.1) -
          POUMM::likPOUMMGivenTreeVTips(
            z = traits.a.1[1,],
            tree = tree.a,
            alpha = model.a.1$H[1,1,1],
            theta = model.a.1$Theta[1,1],
            sigma = model.a.1$Sigma_x[1,1,1],
            sigmae = model.a.1$Sigmae_x[1,1,1],
            g0 = a.X0[1])
        ) < EPS)
  expect_true(
    abs(PCMLik(traits.a.2, tree.a,
        model.a.2) -
          POUMM::likPOUMMGivenTreeVTips(
            z = traits.a.2[1,],
            tree = tree.a,
            alpha = model.a.2$H[1,1,1],
            theta = model.a.2$Theta[1,1],
            sigma = model.a.2$Sigma_x[1,1,1],
            sigmae = model.a.2$Sigmae_x[1,1,1],
            g0 = a.X0[2])
        ) < EPS)
  expect_true(
    abs(PCMLik(traits.a.3, tree.a,
        model.a.3) -
          POUMM::likPOUMMGivenTreeVTips(
            z = traits.a.3[1,],
            tree = tree.a,
            alpha = model.a.3$H[1,1,1],
            theta = model.a.3$Theta[1,1],
            sigma = model.a.3$Sigma_x[1,1,1],
            sigmae = model.a.3$Sigmae_x[1,1,1],
            g0 = a.X0[3])) < EPS)
})

context(ctx <- "OU: R=1/k=1/N=400")

# number of tips
N <- 400

# tree with one regime
tree.a <- rtree(N) # pbtree(n=N, scale=1)
PCMTreeSetDefaultRegime(tree.a, model.a.1)

#tree.a$edge.regime <- rep("a", length(tree.a$edge.length))

# generate traits
traits.a.1 <- PCMSim(tree.a, model.a.1, 0, verbose=TRUE)
traits.a.2 <- PCMSim(tree.a, model.a.2, 0, verbose=TRUE)
traits.a.3 <- PCMSim(tree.a, model.a.2, 0, verbose=TRUE)

test_that(paste(ctx, "Match univariate likelihood from patherit regime a"), {
  expect_true(
    abs(PCMLik(traits.a.1, tree.a,
              model.a.1) -
          POUMM::likPOUMMGivenTreeVTips(
            z = traits.a.1[1,],
            tree = tree.a,
            alpha = model.a.1$H[1,1,1],
            theta = model.a.1$Theta[1,1],
            sigma = model.a.1$Sigma_x[1,1,1],
            sigmae = model.a.1$Sigmae_x[1,1,1],
            g0 = a.X0[1])
    ) < EPS)
  expect_true(
    abs(PCMLik(traits.a.2, tree.a,
              model.a.2) -
          POUMM::likPOUMMGivenTreeVTips(
            z = traits.a.2[1,],
            tree = tree.a,
            alpha = model.a.2$H[1,1,1],
            theta = model.a.2$Theta[1,1],
            sigma = model.a.2$Sigma_x[1,1,1],
            sigmae = model.a.2$Sigmae_x[1,1,1],
            g0 = a.X0[2])
    ) < EPS)
  expect_true(
    abs(PCMLik(traits.a.3, tree.a,
              model.a.3) -
          POUMM::likPOUMMGivenTreeVTips(
            z = traits.a.3[1,],
            tree = tree.a,
            alpha = model.a.3$H[1,1,1],
            theta = model.a.3$Theta[1,1],
            sigma = model.a.3$Sigma_x[1,1,1],
            sigmae = model.a.3$Sigmae_x[1,1,1],
            g0 = a.X0[3])) < EPS)
})

context(ctx <- "OU: R=1/k=3/N=2")

traits.a.123 <- PCMSim(tree.a, model.a.123, c(0,0,0), verbose=TRUE)

test_that(paste(ctx, "Match multivariate likelihood of independent traits regime a"), {
  expect_true(
    abs(PCMLik(traits.a.123, tree.a, model.a.123) -
          (POUMM::likPOUMMGivenTreeVTips(
            traits.a.123[1,],
            tree.a,
            model.a.123$H[1,1,1],
            model.a.123$Theta[1,1],
            model.a.123$Sigma_x[1,1,1],
            model.a.123$Sigmae_x[1,1,1],
            a.X0[1]) +

             POUMM::likPOUMMGivenTreeVTips(
               traits.a.123[2,],
               tree.a,
               model.a.123$H[2,2,1],
               model.a.123$Theta[2,1],
               model.a.123$Sigma_x[2,2,1],
               model.a.123$Sigmae_x[2,2,1],
               a.X0[2]) +

             POUMM::likPOUMMGivenTreeVTips(traits.a.123[3,],
                              tree.a,
                              model.a.123$H[3,3,1],
                              model.a.123$Theta[3,1],
                              model.a.123$Sigma_x[3,3,1],
                              model.a.123$Sigmae_x[3,3,1],
                              a.X0[3]))
        ) < EPS)
})



context(ctx <- "OU: R=1/k=3/N=400")

traits.a.123 <- PCMSim(tree.a, model.a.123, c(0,0,0), verbose=TRUE)

OUlik <- PCMLik(traits.a.123, tree.a, model.a.123)
POUMMlik <- (POUMM::likPOUMMGivenTreeVTips(
  traits.a.123[1,],
  tree.a,
  model.a.123$H[1,1,1],
  model.a.123$Theta[1,1],
  model.a.123$Sigma_x[1,1,1],
  model.a.123$Sigmae_x[1,1,1],
  a.X0[1]) +
    POUMM::likPOUMMGivenTreeVTips(
      traits.a.123[2,],
      tree.a,
      model.a.123$H[2,2,1],
      model.a.123$Theta[2,1],
      model.a.123$Sigma_x[2,2,1],
      model.a.123$Sigmae_x[2,2,1],
      a.X0[2]) +
    POUMM::likPOUMMGivenTreeVTips(
      traits.a.123[3,],
      tree.a,
      model.a.123$H[3,3,1],
      model.a.123$Theta[3,1],
      model.a.123$Sigma_x[3,3,1],
      model.a.123$Sigmae_x[3,3,1],
      a.X0[3]))


cat('OU likelihood=',OUlik,'\n')
cat('POUMM likelihood=',POUMMlik,'\n')

## Calculate likelihood
test_that(paste(ctx, "Match multivariate likelihood of independent traits regime a"), {
  expect_true(abs(OUlik - POUMMlik) < EPS)
})

PCMTreeSetDefaultRegime(tree.a, model.b.123)
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

  PCMTreeSetDefaultRegime(tree.a, model.a.123)
  test_that("a.123",
            expect_equal(PCMLik(traits.a.123, tree.a, model.a.123),
                         PCMLik(traits.a.123, tree.a, model.a.123,
                               metaI = PCMInfoCpp(X = traits.a.123[, 1:length(tree.a$tip.label)],
                                                     tree = tree.a,
                                                     model.a.123))))

  test_that("ab.123",
            expect_equal(PCMLik(traits.ab.123, tree.ab.singles, model.ab.123),
                         PCMLik(traits.ab.123, tree.ab.singles, model.ab.123,
                               metaI = PCMInfoCpp(X = traits.ab.123[, 1:length(tree.ab.singles$tip.label)],
                                                     tree = tree.ab.singles,
                                                     model.ab.123))))


  values <- traits.ab.123[, 1:length(tree.ab.singles$tip.label)]
  values[sample(x=1:length(values), 88)] <- NA

  metaIR <- PCMInfo(X = values,
                     tree = tree.ab.singles,
                     model.ab.123)
  metaICpp <- PCMInfoCpp(X = values,
                               tree = tree.ab.singles,
                               model.ab.123)
  test_that("ab.123 with missing values",
            expect_equal(PCMLik(values, tree.ab.singles, model.ab.123,
                               metaI = metaIR),
                         PCMLik(values, tree.ab.singles, model.ab.123,
                               metaI = metaICpp)))

  print(PCMLik(traits.ab.123, tree.ab.singles, model.ab.123,
               metaI = PCMInfoCpp(
                 X = traits.ab.123[, 1:length(tree.ab.singles$tip.label)],
                 tree = tree.ab.singles,
                 model.ab.123)))

  print(PCMLik(values, tree.ab.singles, model.ab.123,
              metaI = PCMInfoCpp(X = values,
                                    tree = tree.ab.singles,
                                    model.ab.123)))


  if(require(microbenchmark)) {
    cat("microbenchmark test")

    options(PCMBase.PCMLmr.mode=11)
    print(microbenchmark(
      PCMLik(values, tree.ab.singles, model.ab.123, metaI = metaIR),
      PCMLik(values, tree.ab.singles, model.ab.123, metaI = metaICpp),
      times = 10
    ))

    options(PCMBase.PCMLmr.mode=21)
    print(microbenchmark(
      PCMLik(values, tree.ab.singles, model.ab.123, metaI = metaICpp)
    ))
  }
}

if(require(OUwie)) {
  data(tworegime)

  #Calculate the likelihood based on known values of
  #alpha, sigma^2, and theta:
  alpha=c(0.5632459,0.1726052)
  sigma.sq=c(0.1064417,0.3461386)
  theta=c(1.678196,0.4185894)

  OUwie.fixed(tree,trait,model=c("OUMVA"), simmap.tree=FALSE, scaleHeight=FALSE,
              clade=NULL, alpha=alpha,sigma.sq=sigma.sq,theta=theta)


}

