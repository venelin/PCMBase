library(testthat)
context("PCMParam")

library(PCMBase)

if(PCMBaseIsADevRelease(numVersionComponents = 3)) {

  load("testobjects.RData")

  set.seed(1)

  k <- PCMNumTraits(model.ab.123)
  R <- PCMNumRegimes(model.ab.123)
  randVecs1 <- PCMParamRandomVecParams(o = model.ab.123, k = k, R = R, n = 10)

  randVecs2 <- PCMParamRandomVecParams(
    o = model.ab.123, k = k, R = R, n = 10,
    X = traits.ab.123[, seq_len(PCMTreeNumTips(tree.ab))], nUseData = 5)

  randVecs3 <- PCMParamRandomVecParams(
    o = model.ab.123, k = k, R = R, n = 10,
    X = traits.ab.123[, seq_len(PCMTreeNumTips(tree.ab))],
    tree = tree.ab,
    nUseData = 5)

  test_that(
    "First 5 random generated vectors have the grand-mean set for X0",
    expect_equal(
      randVecs2[1:5, 1:3],
      matrix(
        rowMeans(traits.ab.123[, seq_len(PCMTreeNumTips(tree.ab))],
                 na.rm = TRUE),
        nrow = 5, ncol = 3, byrow = TRUE)))

}

