library(ape)
library(testthat)
library(PCMBase)
library(abind)

set.seed(1)

MVTlik <- function(values, tree, model) {
  mean <- model$X0
  sigma <- model$Sigmae[,,1]
  sum(apply(values, 2,
            function(xi) mvtnorm::dmvnorm(x = xi, mean = mean,
                                          sigma = sigma, log = TRUE )))
}


#numerical zero:
EPS <- 10^-8

# number of regimes
R <- 2
# number of traits
k <- 3

# rate mtrix of transition from one regime to another
Q <- matrix(c(-1, 1, 1, -1), R, R)
colnames(Q) <- rownames(Q) <- letters[1:R]

a.X0 <- c(12, 4, 3)
a.Sigmae2 <- rbind(c(.2, 0, 0),
                   c(0, .3, 0),
                   c(0, 0, .4))

# First, specify and sigmae2 parameters for each regime.
# Then we use the abind function to stack the parameters into arrays which's first
# dimension is the regime

Sigmae <- abind(a.Sigmae2, along=3, new.names=list(x=NULL, y=NULL, regime=c('a')))

# regime 'a', traits 1, 2 and 3

model.a.123 <- PCM("White", k = 3, regimes = "a",
                   params = list(X0 = a.X0,
                                 Sigmae = Sigmae[,, "a", drop = FALSE]))


################ 1st Validation ######################################################

context(ctx <- "R=1/k=3/N=2")

# number of tips
N <- 2

# tree with one regime
tree.a <- rtree(N) # phytools::pbtree(n=N, scale=1)
PCMSetDefaultRegime(tree.a, model.a.123)

# generate traits

traits.a.123 <- PCMSim(tree.a, model.a.123, model.a.123$X0, verbose=TRUE)

whitelik = PCMLik(X = traits.a.123[,1:N],
      tree = tree.a,
      model = model.a.123)

mvtlik <- MVTlik(traits.a.123[,1:N],
                tree = tree.a,
                model = model.a.123)

cat('White likelihood=',whitelik,'\n')
cat('mvtnorm likelihood=',mvtlik,'\n')

test_that(paste(ctx, "Match multivariate likelihood of independent traits regime a"), {
  expect_true(abs(whitelik - mvtlik) < EPS)
})

############################# Change Number of tips

context(ctx <- "R=1/k=3/N=400")

# number of tips
N <- 400

tree.a <- rtree(N) # phytools::pbtree(n=N, scale=1)
PCMSetDefaultRegime(tree.a, model.a.123)
# generate traits

traits.a.123 <- PCMSim(tree.a, model.a.123, model.a.123$X0, verbose=TRUE)

# tree with one regime

whitelik = PCMLik(X = traits.a.123,
              tree = tree.a,
              model = model.a.123)

whitelik2 <- PCMLik(X = traits.a.123,
                    tree = tree.a,
                    model = PCM("BM", k = 3, regimes = "a", params = model.a.123))
mvtlik = MVTlik(traits.a.123[, 1:N],
                tree = tree.a,
                model = model.a.123)

cat('white likelihood=',whitelik,'\n')
cat('white likelihood 2=',whitelik2,'\n')
cat('mvtlik likelihood=',mvtlik,'\n')

test_that(paste(ctx, "Match multivariate likelihood of independent traits regime a"), {
  expect_true(abs(whitelik - mvtlik) < EPS)
})

if(require(PCMBaseCpp)) {
  cat("Testing PCMBaseCpp on White:\n")

  whitelik <- PCMLik(tree = tree.a, model = model.a.123,
                      metaI = PCMInfoCpp(X = traits.a.123[, 1:N], tree = tree.a, model = model.a.123))

  test_that("Match white likelihoods", expect_true(abs(whitelik-mvtlik) < EPS))

  values <- traits.a.123[, 1:length(tree.a$tip.label)]
  values[sample(x=1:length(values), 20)] <- NA

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

