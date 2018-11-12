library(testthat)
context("PCMLik, OU")

library(PCMBase)

if(PCMBaseIsADevRelease(numVersionComponents = 3)) {

  library(mvtnorm)

  load("testobjects.RData")



  # test likelihood
  test_that(
    "Single trait log-likelihood, regime a", {
      expect_equivalent(PCMLik(traits.a.1, tree.a, model.a.1), -91.015331180479)
      expect_equivalent(PCMLik(traits.a.2, tree.a, model.a.2), -60.0600001079255)
      expect_equivalent(PCMLik(traits.a.3, tree.a, model.a.3), -527.311935254892)
    })


  test_that(
    "Triple-trait log-likelihood of independent traits, regime a", {
      expect_equivalent(PCMLik(traits.a.123, tree.a, model.a.123), -205.993838713138)
    })


  MeanVec <- PCMMean(tree.a, model.a.123, model.a.123$X0)
  VarMat <- PCMVar(tree.a, model.a.123)

  test_that(
    "Triple-trait log-likelihood, of independent traits, regime a, using PCMMean and PCMVar", {
      expect_equivalent(
        dmvnorm(as.vector(traits.a.123[, 1:PCMTreeNumTips(tree.a)]),
                as.vector(MeanVec), VarMat, log = TRUE),
        -205.993838713138
      )
    }
  )
}
