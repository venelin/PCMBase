.RunPCMBaseTests <- Sys.getenv("RunPCMBaseTests") == "yes"

if(.RunPCMBaseTests) {

# Test for ScalarDiagonal_Sigma_x matrices - this causes a problem with the default
# arma::eig_gen implementation.
#
library(ape)
library(testthat)
library(PCMBase)
library(PCMFit)
library(abind)
library(data.table)
library(phytools)

set.seed(2)

# number of regimes
R <- 3
# number of traits
k <- 2

# number of tips
N <- 200

tree.a <- pbtree(n=N, scale=1, b = 1, d = 0.4)
PCMTreeSetRegimes(tree.a, c(322, 513), regimes = c("a", "b", "c"))
PCMTreeSetLabels(tree.a)

tree.a$edge.jump <- rep(0, nrow(tree.a$edge))

## Try to reproduce the error outside MixedGaussian
listParameterizationsBM2 <- list(
  X0 = list(c("VectorParameter", "_Global")),
  Sigma_x = list(c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
                 c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
                 c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal")),

  Sigmae_x = list(c("MatrixParameter", "_Zeros"))
)

listParameterizationsOU2 <- list(
  X0 = list(c("VectorParameter", "_Global")),
  H = list(c("MatrixParameter", "_Zeros")),
  Theta = list(c("VectorParameter", "_Zeros")),
  Sigma_x = list(c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
                 c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
                 c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal")),

  Sigmae_x = list(c("MatrixParameter", "_Zeros"))
)


PCMGenerateParameterizations(structure(0.0, class="BM"), listParameterizations = listParameterizationsBM2)
PCMGenerateParameterizations(structure(0.0, class="OU"), listParameterizations = listParameterizationsOU2)

model <- PCM("BM__Global_X0__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Zeros_Sigmae_x",
              modelTypes = c("BM__Global_X0__ScalarDiagonal_WithNonNegativeDiagonal_Sigma_x__Zeros_Sigmae_x",
                             "BM__Global_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__Zeros_Sigmae_x",
                             "BM__Global_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Zeros_Sigmae_x"),
              k = k,
              regimes = letters[1:3])

vecModelRandom <- round(PCMParamRandomVecParams(model), 1)
modelRandom <- model
PCMParamLoadOrStore(modelRandom, vecModelRandom, offset = 0, load=TRUE)
model <- modelRandom

traits <- PCMSim(tree.a, model, X0 = model$X0, verbose=TRUE)

options(PCMBase.Value.NA = -1e20)
values <- traits[, 1:length(tree.a$tip.label)]

likR <- PCMLik(values, tree.a, model)

if(require(PCMBaseCpp)) {
  metaInfo <- PCMInfoCpp(values, tree.a, model)
  likFun2 <- PCMCreateLikelihood(values, tree.a, model, metaInfo)

  print(likFun2(p = PCMParamGetShortVector(model)))
  print(likR)
  test_that("Equal R and PCMBaseCpp values with scalar diagonal Sigma matrix",
            expect_equal(likFun2(p = PCMParamGetShortVector(model)), likR))

}



}
