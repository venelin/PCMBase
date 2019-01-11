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
    o = model.ab.123, k = k, R = R, n = 10)

  randVecs3 <- PCMParamRandomVecParams(
    o = model.ab.123, k = k, R = R, n = 10)

}

