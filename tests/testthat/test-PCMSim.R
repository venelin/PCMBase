library(testthat)
context("PCMSim")

library(PCMBase)

if(PCMBaseIsADevRelease(numVersionComponents = 3)) {


  load("testobjects.RData")

  set.seed(1)

  # generate traits
  traits.a.1 <- PCMSim(tree.a, model.a.1, 0, verbose=TRUE)

  test_that(
    "Simulated single-trait data on a single regime tree...", {
      expect_identical(dim(traits.a.1), c(1L, PCMTreeNumNodes(tree.a)))
    }
  )

  traits.a.2 <- PCMSim(tree.a, model.a.2, 0, verbose=TRUE)
  traits.a.3 <- PCMSim(tree.a, model.a.2, 0, verbose=TRUE)

  traits.a.123 <- PCMSim(tree.a, model.a.123, c(0,0,0), verbose=TRUE)

  test_that(
    "Simulated triple-trait data on a single regime tree...", {
      expect_identical(dim(traits.a.123), c(3L, PCMTreeNumNodes(tree.a)))
    }
  )

  traits.ab.123 <- PCMSim(tree.ab, model.ab.123, c(0,0,0), verbose=TRUE)
  test_that(
    "Simulated triple-trait data on a two-regime tree...", {
      expect_identical(dim(traits.ab.123), c(3L, PCMTreeNumNodes(tree.ab)))
    }
  )

  # save testobjects:
  #
  # save(a.H, a.Sigma_x, a.Sigmae_x, b.H, b.Sigma_x, b.Sigmae_x, model_MixedGaussian_ab, model_MixedGaussian_ab_globalSigmae_x, model.a.1, model.a.1.Omitted_X0, model.a.123, model.a.123.Omitted_X0, model.a.123.Omitted_X0__bSigmae_x, model.a.123.Omitted_X0__Omitted_Sigmae_x, model.a.2, model.a.3, model.ab.123, model.ab.123.bSigmae_x, model.b.123, model.b.123.Omitted_X0, model.b.123.Omitted_X0__Omitted_Sigmae_x, tree.a, tree.ab, Theta, a.Theta, a.X0, b.Theta, b.X0, H, N, Sigma_x, Sigmae_x, traits.a.1, traits.a.123, traits.a.2, traits.a.3, traits.ab.123, model.a.1, file = "tests/testthat/testobjects.RData")

}
