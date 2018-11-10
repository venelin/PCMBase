.RunPCMBaseTests <- Sys.getenv("RunPCMBaseTests") == "yes"

if(.RunPCMBaseTests) {

  library(ape)
  library(testthat)
  library(abind)
  library(PCMBase)

  set.seed(1)

  MVTlik <- function(values, tree, model) {
    mean <- model$X0
    sigma <- as.matrix(model$Sigmae_x[,,1]) %*% t(as.matrix(model$Sigmae_x[,,1]))
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
  a.Sigmae_x <- rbind(c(.2, 0, 0),
                      c(0, .3, 0),
                      c(0, 0, .4))

  # First, specify and sigmae2 parameters for each regime.
  # Then we use the abind function to stack the parameters into arrays which's first
  # dimension is the regime

  Sigmae_x <- abind(a.Sigmae_x, along=3, new.names=list(x=NULL, y=NULL, regime=c('a')))

  # regime 'a', traits 1, 2 and 3

  model.a.123 <- PCM("White", k = 3, regimes = "a",
                     params = list(X0 = a.X0,
                                   Sigmae_x = Sigmae_x[,, "a", drop = FALSE]))


  ################ 1st Validation ######################################################

  context(ctx <- "R=1/k=3/N=2")

  # number of tips
  N <- 2

  # tree with one regime
  tree.a <- rtree(N) # phytools::pbtree(n=N, scale=1)
  PCMTreeSetDefaultRegime(tree.a, model.a.123)

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
  PCMTreeSetDefaultRegime(tree.a, model.a.123)
  # generate traits

  traits.a.123 <- PCMSim(tree.a, model.a.123, model.a.123$X0, verbose=TRUE)

  # tree with one regime

  whitelik = PCMLik(X = traits.a.123,
                    tree = tree.a,
                    model = model.a.123)

  mvtlik = MVTlik(traits.a.123[, 1:N],
                  tree = tree.a,
                  model = model.a.123)

  cat('white likelihood=',whitelik,'\n')
  cat('mvtlik likelihood=',mvtlik,'\n')

  test_that(paste(ctx, "Match multivariate likelihood of independent traits regime a"), {
    expect_true(abs(whitelik - mvtlik) < EPS)
  })



}
