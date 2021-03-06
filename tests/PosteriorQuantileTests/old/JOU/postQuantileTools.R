#' function that converts a vector of JOU parameters to a list of matrices
#' @param a vector of OU parameters of length 14
#' @return a list of matrices.

fromVector <- function(l){

  Alpha = diag(x=c(l[1],l[2]), 2, 2)
  Theta = as.vector(c(l[3],l[4]))
  Sigma = matrix(c(l[5],l[6],l[6],l[7]), 2, 2,byrow = TRUE)
  Sigmae = diag(x=c(l[8],l[9]), 2, 2)
  mj = as.vector(c(l[10],l[11]))
  Sigmaj = matrix(c(l[12],l[13],l[13],l[14]), 2, 2,byrow = TRUE)

  param.m = list('Alpha' = Alpha,'Theta'= Theta,'Sigma' = Sigma,'Sigmae'=Sigmae,'mj'=mj,'Sigmaj'=Sigmaj)

}

#' function that generates the seeds for a multivariate JOU porcess according to specified parameter distributions.
#' @param the number of replications
#' @return a list of matrices.

generate.seeds <- function(L){

  dAlpha = lapply(1:2,function(i) rexp(L,1))
  Alpha.replicates = lapply(1:L, function(i) diag(x=c(dAlpha[[1]][i],dAlpha[[2]][i]), 2, 2))

  Theta1 = rnorm(L,0,1)
  Theta2 = rnorm(L,2,1)
  Theta.replicates = lapply(1:L, function(i) c(Theta1[i],Theta2[i]))

  mj1 = rnorm(L,0.2,0.1)
  mj2 = rnorm(L,-0.2,0.1)
  mj.replicates = lapply(1:L, function(i) c(mj1[i],mj2[i]))

  Sigma.replicates = lapply(1:L, function(i) {
    determ <- -Inf
    while(determ<0) {
      dSigma = sapply(1:2,function(i) rexp(1,1))
      offdSigma = rnorm(1,0.2,0.1)

      mat <- matrix(c(dSigma[1],offdSigma,offdSigma,dSigma[2]), 2, 2,
                    byrow = TRUE)
      determ <- det(mat)
    }
    mat
  })


  Sigmaj.replicates = lapply(1:L, function(i) {
    determ <- -Inf
    while(determ<0) {
      dSigmaj = sapply(1:2,function(i) rexp(1,5))
      offdSigmaj = rnorm(1,0.04,0.025)

      mat <- matrix(c(dSigmaj[1],offdSigmaj,offdSigmaj,dSigmaj[2]), 2, 2,
                    byrow = TRUE)
      determ <- det(mat)
    }
    mat
  })

  dSigmae = lapply(1:2,function(i) rexp(L,10))
  Sigmae.replicates = lapply(1:L, function(i) diag(x=c(dSigmae[[1]][i],dSigmae[[2]][i]), 2, 2))

  return(list('Alpha' = Alpha.replicates,'Theta'=Theta.replicates,'Sigma'=Sigma.replicates,
              'Sigmae'=Sigmae.replicates,'mj'=mj.replicates,'Sigmaj'=Sigmaj.replicates))

}

#' function that generates a JOU model
#' @param parameters is a list of matrices containing the parameters Alpha, Theta, Sigma, Sigmae, mj and Sigmaj
#' @param binary vector of length equal to the number of edges in the tree indicating a jump or not at each edge
#' @return JOU model.

generate.model <- function(parameters,xi){

  Alpha1 = parameters[['Alpha']]
  Theta1 = parameters[['Theta']]
  Sigma1 = parameters[['Sigma']]
  Sigmae1 = parameters[['Sigmae']]
  mj1 = parameters[['mj']]
  Sigmaj1 = parameters[['Sigmaj']]


  Alpha <- abind::abind(Alpha1, Alpha1, along=-1, new.names=list(regime=c('a','a2'), x=NULL, y=NULL))
  Theta <- abind::abind(Theta1, Theta1, along=-1, new.names=list(regime=c('a', 'a2'), xy=NULL))
  Sigma <- abind::abind(Sigma1, Sigma1, along=-1, new.names=list(regime=c('a','a2'), x=NULL, y=NULL))
  Sigmae <- abind::abind(Sigmae1, Sigmae1, along=-1, new.names=list(regime=c('a','a2'), x=NULL, y=NULL))
  mj <- abind::abind(mj1, mj1, along=-1, new.names=list(regime=c('a', 'a2'), xy=NULL))
  Sigmaj <- abind::abind(Sigmaj1, Sigmaj1, along=-1, new.names=list(regime=c('a','a2'), x=NULL, y=NULL))


  model <- list(X0 = c(0,0),
                Alpha=Alpha['a',,,drop=FALSE],
                Theta=Theta['a',,drop=FALSE],
                Sigma=Sigma['a',,,drop=FALSE],
                Sigmae=Sigmae['a',,,drop=FALSE],
                mj=mj['a',,drop=FALSE],
                Sigmaj=Sigmaj['a',,,drop=FALSE],
                xi=xi)
  class(model) <- 'JOU'

  return(model)

}

#' function that generates the traits on a tree
#' @param tree
#' @param JOU model
#' @return a list containing the tree and the traits.

generate.data <-function(tree,model){

  traits <- mvsim(tree, model, c(0,0), verbose=TRUE)

  return(list('tree'=tree,'traits'=traits))
}

#' function that computes the prior given the distribution of the parameters
#' @param a list of matrices containing the parameters of the JOU model: Alpha, Theta, Sigma, Sigmae, mj and Sigmaj
#' @return a single value for the prior or -Inf if Sigma or Sigmaj matrix is not positive semi-definite

compute.prior <- function(parameters){

  if (det(parameters[['Sigma']])<0 | det(parameters[['Sigmaj']])<0){
    return (-Inf)
  }
  return (sum(c(sum(dexp(parameters[['Alpha']],log=TRUE)),sum(dnorm(parameters[['Theta']],c(0,2),c(1,1),log=TRUE)),
                sum(dexp(diag(parameters[['Sigma']]),log=TRUE)),dnorm(parameters[['Sigma']][1,2],0.2,0.1,log=TRUE),
                sum(dexp(parameters[['Sigmae']],log=TRUE)),sum(dnorm(parameters[['mj']],c(0.2,-0.2),c(0.1,0.1),log=TRUE)),
                sum(dexp(diag(parameters[['Sigmaj']]),rate=5,log=TRUE)),dnorm(parameters[['Sigmaj']][1,2],0.04,0.025,log=TRUE))))

}

#' function that computes a value proportional to the log probability density from which the MCMC can sample.
#' @param vector of parameters in the JOU model
#' @param list containing the tree and the traits
#' @return the sum of the likelihood of the JOU parameters and the prior value.

posterior <- function(parameters,extra){

  tree = extra[[1]]
  traits = extra[[2]]
  xi = extra[[3]]
  parameters.m = fromVector(parameters)
  model = generate.model(parameters.m,xi)
  prior_value = compute.prior(parameters.m)

  if (is.infinite(prior_value)){
    sum = -Inf
  }else{
    lik_value = mvlik(X = traits$values + traits$errors, tree = tree, model = model)
    sum = prior_value + lik_value[1]
  }
  return (sum)

}

#' function that generates the MCMC samples.
#' @param posterior function
#' @param number of samples
#' @param the tree
#' @param the traits of the tree
#' @return a list of components returned by the built-in R MCMC function.

generate.samples <- function(p,n,tree,traits,xi){

  samples = MCMC(p, n,init = c(1,1,0,2,1,0.2,1,0.1,0.1,0.2,-0.2,5,0.04,5), scale = rep(1,14),adapt = TRUE, acc.rate = 0.234,
                 list = TRUE, extra = list(tree,traits,xi))

  return (samples)
}

#' function that computes the posterior quantiles for all JOU parameters.
#' @param the samples generated from the MCMC in the form of coda objects
#' @param a list of the seeds
#' @return a vector of posterior quantiles with length equal to the number of JOU parameters.

compute.postquant <- function(samples,initial){

  num.param = ncol(samples)
  num.samples = nrow(samples)

  res = sapply(1:num.param, function(i) (sum(samples[,i]<initial[i])) / num.samples)
  return(res)
}
