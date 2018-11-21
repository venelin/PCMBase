library(abind)
library(PCMBase)

# Generate some parameters and models:

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


# regime 'a', trait 1
model.a.1 <- PCM("OU", k = 1, regimes = "a",
                 params = list(
                   X0 = a.X0[1],
                   H = H[1,1,'a',drop=FALSE],
                   Theta = Theta[1,'a',drop=FALSE],
                   Sigma_x = Sigma_x[1,1,'a',drop=FALSE],
                   Sigmae_x = Sigmae_x[1,1,'a',drop=FALSE]))

# regime 'a', trait 2
model.a.2 <- PCM("OU", k = 1, regimes = "a",
                 params = list(
                   X0 = a.X0[2],
                   H = H[2,2,'a',drop=FALSE],
                   Theta = Theta[2,'a',drop=FALSE],
                   Sigma_x = Sigma_x[2,2,'a',drop=FALSE],
                   Sigmae_x = Sigmae_x[2,2,'a',drop=FALSE]))

# regime 'a', trait 3
model.a.3 <- PCM("OU", k = 1, regimes = "a",
                 params = list(
                   X0 = a.X0[3],
                   H = H[3,3,'a',drop=FALSE],
                   Theta = Theta[3,'a',drop=FALSE],
                   Sigma_x = Sigma_x[3,3,'a',drop=FALSE],
                   Sigmae_x = Sigmae_x[3,3,'a',drop=FALSE]))

# regime 'a', traits 1, 2 and 3
model.a.123 <- PCM("OU", k = 3, regimes = "a",
                   params = list(
                     X0 = a.X0,
                     H = H[,,'a',drop=FALSE],
                     Theta = Theta[,'a',drop=FALSE],
                     Sigma_x = Sigma_x[,,'a',drop=FALSE],
                     Sigmae_x = Sigmae_x[,,'a',drop=FALSE]))



# regime 'b', traits 1, 2 and 3
model.b.123 <- PCM("OU", k = 3, regimes = "b",
                   params = list(
                     X0 = b.X0,
                     H = H[,,'b',drop=FALSE],
                     Theta = Theta[,'b',drop=FALSE],
                     Sigma_x = Sigma_x[,,'b',drop=FALSE],
                     Sigmae_x = Sigmae_x[,,'b',drop=FALSE]))

# regimes 'a' and 'b', traits 1, 2 and 3
model.ab.123 <- PCM("OU", k = 3, regimes = c("a", "b"),
                    params = list(
                      X0 = a.X0,
                      H = H[,,,drop=FALSE],
                      Theta = Theta[,,drop=FALSE],
                      Sigma_x = Sigma_x[,,,drop=FALSE],
                      Sigmae_x = Sigmae_x[,,,drop=FALSE]))


# generate some  models
tableParametrizationsOU <- PCMTableParameterizations(structure(0.0, class="OU"))
PCMGenerateParameterizations(
      model = structure(0.0, class="OU"),
      tableParameterizations = tableParametrizationsOU[
        sapply(X0, function(type)
          identical(type, c("VectorParameter", "_Global")) ||
            identical(type, c("VectorParameter", "_Omitted"))
        ) &
          sapply(H, function(type)
            identical(type, c("MatrixParameter"))) &
          sapply(Theta, function(type)
            identical(type, "VectorParameter") )
        ])

tableParametrizationsBM <- PCMTableParameterizations(structure(0.0, class="BM"))
PCMGenerateParameterizations(
      model = structure(0.0, class="BM"),
      tableParameterizations = tableParametrizationsBM[
        sapply(X0, function(type)
          identical(type, c("VectorParameter", "_Global")) ||
            identical(type, c("VectorParameter", "_Omitted")) ), ])


# regime 'a', trait 1
model.a.1.Omitted_X0 <- PCM(
  "OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x",
  k = 1,
  regimes = "a",
  params = list(H=H[1,1,'a',drop=FALSE],
                Theta=Theta[1,'a',drop=FALSE],
                Sigma_x=Sigma_x[1,1,'a',drop=FALSE],
                Sigmae_x=Sigmae_x[1,1,'a',drop=FALSE]))


# regime 'a', traits 1, 2 and 3
model.a.123.Omitted_X0 <- PCM(
  "OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x",
  k = 3,
  regimes = "a",
  params = list(H=H[,,'a',drop=FALSE],
                Theta=Theta[,'a',drop=FALSE],
                Sigma_x=Sigma_x[,,'a',drop=FALSE],
                Sigmae_x=Sigmae_x[,,'a',drop=FALSE]))


# regime 'b', traits 1, 2 and 3
model.b.123.Omitted_X0 <- PCM(
  "OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x",
  k = 3L,
  regimes = "b",
  params = list(H=H[,,'b',drop=FALSE],
                Theta=Theta[,'b',drop=FALSE],
                Sigma_x=Sigma_x[,,'b',drop=FALSE],
                Sigmae_x=Sigmae_x[,,'b',drop=FALSE]))

model_MixedGaussian_ab <- MixedGaussian(
  k = 3L,
  modelTypes = "OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x",
  mapping = c(a = 1, b = 1),
  className = "MixedGaussian_ab",
  Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted")))

PCMParamSetByName(model_MixedGaussian_ab,
                  list(X0 = a.X0,
                       a = model.a.123.Omitted_X0,
                       b = model.b.123.Omitted_X0))


# regime 'a', traits 1, 2 and 3
model.a.123.Omitted_X0__bSigmae_x <- PCM(
  "OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x", k = 3, regimes = "a",
  params = list(H=H[,,'a',drop=FALSE],
                Theta=Theta[,'a',drop=FALSE],
                Sigma_x=Sigma_x[,,'a',drop=FALSE],
                Sigmae_x=Sigmae_x[,,'b',drop=FALSE]))


# regimes 'a' and 'b', traits 1, 2 and 3
model.ab.123.bSigmae_x <- PCM(
  "OU", k = 3, regimes = c("a", "b"),
  params = list(X0 = a.X0,
                H=H[,,,drop=FALSE],
                Theta=Theta[,,drop=FALSE],
                Sigma_x=Sigma_x[,,,drop=FALSE],
                Sigmae_x=Sigmae_x[,,c('b', 'b'),drop=FALSE]))

model.a.123.Omitted_X0__Omitted_Sigmae_x <- PCM(
  "OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x", k = 3, regimes = "a",
  params = list(H=H[,,'a',drop=FALSE],
                Theta=Theta[,'a',drop=FALSE],
                Sigma_x=Sigma_x[,,'a',drop=FALSE]))

model.b.123.Omitted_X0__Omitted_Sigmae_x <- PCM(
  "OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x", k = 3, regimes = "b",
  params = list(H=H[,,'b',drop=FALSE],
                Theta=Theta[,'b',drop=FALSE],
                Sigma_x=Sigma_x[,,'b',drop=FALSE]))


model_MixedGaussian_ab_globalSigmae_x <- MixedGaussian(
  k = 3,
  modelTypes = "OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
  mapping = c(a = 1, b = 1),
  className = "MixedGaussian_ab_globalSigmae_x")

PCMParamSetByName(model_MixedGaussian_ab_globalSigmae_x,
                  list(X0 = a.X0,
                       a = model.a.123.Omitted_X0__Omitted_Sigmae_x,
                       b = model.b.123.Omitted_X0__Omitted_Sigmae_x,
                       Sigmae_x=Sigmae_x[,,'b',drop=TRUE]))

set.seed(1)

model.ab.123.MG <- MixedGaussian(
    k = 3,
    modelTypes = c("BM__Omitted_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
                   "OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"),
    mapping = c(a=2L, b=1L), className = "MG")

PCMParamLoadOrStore(model.ab.123.MG,
                    PCMParamRandomVecParams(model.ab.123.MG),
                    offset = 0, k = 3, load = TRUE)


# generate trees:
library(ape)

set.seed(1)
# number of tips
N <- 40

# tree with one regime
tree.a <- rtree(N)

PCMTreeSetDefaultRegime(tree.a, model.a.1)
PCMTreeSetLabels(tree.a)

#PCMTreePlot(tree.a)

tree.ab <- tree.a
PCMTreeSetRegimes(tree.ab, nodes = N + 31, regimes = c("a", "b"))
#PCMTreePlot(tree.ab)


# generate trait values

set.seed(1)

# generate traits
traits.a.1 <- PCMSim(tree.a, model.a.1, 0, verbose=TRUE)

traits.a.2 <- PCMSim(tree.a, model.a.2, 0, verbose=TRUE)
traits.a.3 <- PCMSim(tree.a, model.a.2, 0, verbose=TRUE)

traits.a.123 <- PCMSim(tree.a, model.a.123, c(0,0,0), verbose=TRUE)

traits.ab.123 <- PCMSim(tree.ab, model.ab.123, c(0,0,0), verbose=TRUE)


# save testobjects:
#

if(TRUE) {
  save(
    a.H,
    a.Sigma_x,
    a.Sigmae_x,
    b.H,
    b.Sigma_x,
    b.Sigmae_x,
    model_MixedGaussian_ab,
    model_MixedGaussian_ab_globalSigmae_x,
    model.a.1,
    model.a.1.Omitted_X0,
    model.a.123,
    model.a.123.Omitted_X0,
    model.a.123.Omitted_X0__bSigmae_x,
    model.a.123.Omitted_X0__Omitted_Sigmae_x,
    model.a.2,
    model.a.3,
    model.ab.123,
    model.ab.123.bSigmae_x,
    model.b.123,
    model.b.123.Omitted_X0,
    model.b.123.Omitted_X0__Omitted_Sigmae_x,
    tree.a,
    tree.ab,
    Theta,
    a.Theta,
    a.X0,
    b.Theta,
    b.X0,
    H,
    N,
    Sigma_x,
    Sigmae_x,
    traits.a.1,
    traits.a.123,
    traits.a.2,
    traits.a.3,
    traits.ab.123,
    model.a.1,

    file = "../tests/testthat/testobjects.RData")
}
