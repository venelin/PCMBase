fromVector <- function(l){

  Alpha1 = diag(x=c(l[1],l[2]), 2, 2)
  Alpha2 = diag(x=c(l[3],l[4]), 2, 2)
  Theta = as.vector(c(l[5],l[6]))
  Sigma = matrix(c(l[7],l[8],l[8],l[9]), 2, 2,byrow = TRUE)
  Sigmae = diag(x=c(l[10],l[11]), 2, 2)

  param.m = list('Alpha1' = Alpha1,'Alpha2' = Alpha2,'Theta'= Theta,'Sigma' = Sigma,'Sigmae'=Sigmae)

}

# L is the number of replicates
generate.seeds <- function(L){

  dAlpha1 = lapply(1:2,function(i) rexp(L,5))
  Alpha1.replicates = lapply(1:L, function(i) diag(x=c(dAlpha1[[1]][i],dAlpha1[[2]][i]), 2, 2))

  dAlpha2 = lapply(1:2,function(i) rexp(L,1))
  Alpha2.replicates = lapply(1:L, function(i) diag(x=c(dAlpha2[[1]][i],dAlpha2[[2]][i]), 2, 2))

  Theta1 = rnorm(L,0,1)
  Theta2 = rnorm(L,2,1)
  Theta.replicates = lapply(1:L, function(i) c(Theta1[i],Theta2[i]))

  # dSigma = lapply(1:2,function(i) rexp(1,1))
  # offdSigma = rnorm(1,0.2,0.1)

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

  dSigmae = lapply(1:2,function(i) rexp(L,10))
  Sigmae.replicates = lapply(1:L, function(i) diag(x=c(dSigmae[[1]][i],dSigmae[[2]][i]), 2, 2))

  return(list('Alpha1' = Alpha1.replicates,'Alpha2' = Alpha2.replicates,'Theta'=Theta.replicates,
              'Sigma'=Sigma.replicates,'Sigmae'=Sigmae.replicates))

}

generate.model <- function(parameters){

  Alpha1 = parameters[['Alpha1']]
  Alpha2 = parameters[['Alpha2']]
  Theta1 = parameters[['Theta']]
  Sigma1 = parameters[['Sigma']]
  Sigmae1 = parameters[['Sigmae']]


  Alpha1 <- abind::abind(Alpha1, Alpha1, along=-1, new.names=list(regime=c('a','a2'), x=NULL, y=NULL))
  Alpha2 <- abind::abind(Alpha2, Alpha2, along=-1, new.names=list(regime=c('a','a2'), x=NULL, y=NULL))
  Theta <- abind::abind(Theta1, Theta1, along=-1, new.names=list(regime=c('a', 'a2'), xy=NULL))
  Sigma <- abind::abind(Sigma1, Sigma1, along=-1, new.names=list(regime=c('a','a2'), x=NULL, y=NULL))
  Sigmae <- abind::abind(Sigmae1, Sigmae1, along=-1, new.names=list(regime=c('a','a2'), x=NULL, y=NULL))

  model <- list(X0 = c(0,0),
                Alpha1=Alpha1['a',,,drop=FALSE],
                Alpha2=Alpha2['a',,,drop=FALSE],
                Theta=Theta['a',,drop=FALSE],
                Sigma=Sigma['a',,,drop=FALSE],
                Sigmae=Sigmae['a',,,drop=FALSE])
  class(model) <- '2SpOU'

  return(model)

}

generate.data <-function(tree,model){

  traits <- mvsim(tree, model, c(0,0), verbose=TRUE)

  return(list('tree'=tree,'traits'=traits))
}

# parameters are in a matrix form.
compute.prior <- function(parameters){

  if (det(parameters[['Sigma']])<0){
    return (-Inf)
  }
  return (sum(c(sum(dexp(parameters[['Alpha1']],rate=5,log=TRUE)),
                sum(dexp(parameters[['Alpha2']],log=TRUE)),
                sum(dnorm(parameters[['Theta']],c(0,2),c(1,1),log=TRUE)),
                sum(dexp(diag(parameters[['Sigma']]),log=TRUE)),dnorm(parameters[['Sigma']][1,2],0.2,0.1,log=TRUE),
                sum(dexp(parameters[['Sigmae']],log=TRUE)))))

}

# parameters are in a list form
posterior <- function(parameters,extra){

  tree = extra[[1]]
  traits = extra[[2]]
  parameters.m = fromVector(parameters)
  model = generate.model(parameters.m)
  prior_value = compute.prior(parameters.m)

  if (is.infinite(prior_value)){
    sum = -Inf
  }else{
    lik_value = mvlik(X = traits$values + traits$errors, tree = tree, model = model)
    sum = prior_value + lik_value[1]
  }
  return (sum)

}

generate.samples <- function(p,n,tree,traits){

  samples = MCMC(p, n,init = c(5,5,1,1,0,2,1,0.2,1,0.1,0.1), scale = rep(1,11),adapt = TRUE, acc.rate = 0.234,
                 list = TRUE, extra = list(tree,traits))

  return (samples)
}

compute.postquant <- function(samples,initial){

  num.param = ncol(samples)
  num.samples = nrow(samples)

  res = sapply(1:num.param, function(i) (sum(samples[,i]<initial[i])) / num.samples)
  return(res)
}
