# 1st Draft is for OU

################################################MAIN
library(adaptMCMC)

replications = 15
samples_res = list()
post.quant = list()
length(samples_res) = replications
length(post.quant) = replications
param.replicates <- generate.seeds(replications)


for (iter in 1:replications){

  Alpha = param.replicates[['Alpha']][[iter]]   ################## remember to put iter
  Theta = param.replicates[['Theta']][[iter]]
  Sigma = param.replicates[['Sigma']][[iter]]
  Sigmae = param.replicates[['Sigmae']][[iter]]

  param.m = list('Alpha' = Alpha,'Theta'= Theta,'Sigma' = Sigma,'Sigmae'=Sigmae)  ######### or I take out the names

  model = generate.model(param.m)

  nodes = 100
  data = generate.data(nodes,model)
  tree = data$tree
  traits = data$traits

  param.l = c(Alpha[1,1],Alpha[2,2],Theta,Sigma[1,1],Sigma[1,2],Sigma[2,2],Sigmae[1,1],Sigmae[2,2])

  samples_res[[iter]] = generate.samples(posterior,10,tree,traits)

  post.quant[[iter]] = compute.postquant(samples_res[[iter]]$samples,param.l)

}

param.names = c('Alpha11','Alpha22','Theta1','Theta2','Sigma11','Sigma12','Sigma22','Sigmae11','Sigmae22')
num.param = length(param.l)

x = list()
length(x) = num.param

for (j in 1:num.param){

  x[[j]] = sapply(1:replications, function(i) (c(post.quant[[i]][j])))
  filename = paste(param.names[j],"histOU.png")
  png(filename)
  hist(x[[j]], main=paste("Histogram for OU parameter",param.names[j]),xlab = param.names[j]
       ,border="black",col="blue",las=1, breaks=5)
  dev.off()
}

################################################


#input is a list of values
fromVector <- function(l){

  Alpha = diag(x=c(l[1],l[2]), 2, 2)
  Theta = as.vector(c(l[3],l[4]))
  Sigma = matrix(c(l[5],l[6],l[6],l[7]), 2, 2,byrow = TRUE)
  Sigmae = diag(x=c(l[8],l[9]), 2, 2)

  param.m = list('Alpha' = Alpha,'Theta'= Theta,'Sigma' = Sigma,'Sigmae'=Sigmae)

}

# L is the number of replicates
generate.seeds <- function(L){

  dAlpha = lapply(1:2,function(i) rexp(L,1))
  Alpha.replicates = lapply(1:L, function(i) diag(x=c(dAlpha[[1]][i],dAlpha[[2]][i]), 2, 2))

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

  return(list('Alpha' = Alpha.replicates,'Theta'=Theta.replicates,
           'Sigma'=Sigma.replicates,'Sigmae'=Sigmae.replicates))

}

generate.model <- function(parameters){

  Alpha1 = parameters[['Alpha']]
  Theta1 = parameters[['Theta']]
  Sigma1 = parameters[['Sigma']]
  Sigmae1 = parameters[['Sigmae']]


  Alpha <- abind::abind(Alpha1, Alpha1, along=-1, new.names=list(regime=c('a','a2'), x=NULL, y=NULL))
  Theta <- abind::abind(Theta1, Theta1, along=-1, new.names=list(regime=c('a', 'a2'), xy=NULL))
  Sigma <- abind::abind(Sigma1, Sigma1, along=-1, new.names=list(regime=c('a','a2'), x=NULL, y=NULL))
  Sigmae <- abind::abind(Sigmae1, Sigmae1, along=-1, new.names=list(regime=c('a','a2'), x=NULL, y=NULL))

  model <- list(X0 = c(0,0),
                      Alpha=Alpha['a',,,drop=FALSE],
                      Theta=Theta['a',,drop=FALSE],
                      Sigma=Sigma['a',,,drop=FALSE],
                      Sigmae=Sigmae['a',,,drop=FALSE])
  class(model) <- 'OU'

  return(model)

}

generate.data <-function(N,model){

  tree <- phytools::pbtree(n=N, scale=1)
  traits <- mvsim(tree, model, c(0,0), verbose=TRUE)

  return(list('tree'=tree,'traits'=traits))
}

# parameters are in a matrix form.
compute.prior <- function(parameters){

  if (det(parameters[['Sigma']])<0){
    return (-Inf)
  }
  return (sum(c(sum(dexp(parameters[['Alpha']],log=TRUE)),sum(dnorm(parameters[['Theta']],c(0,2),c(1,1),log=TRUE)),
                sum(dexp(diag(parameters[['Sigma']]),log=TRUE)),dnorm(parameters[['Sigma']][1,2],0.2,0.1,log=TRUE),sum(dexp(parameters[['Sigmae']],log=TRUE)))))

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

  samples = MCMC(p, n,init = c(1,1,0,2,1,0.2,1,0.1,0.1), scale = rep(1,9),adapt = TRUE, acc.rate = 0.234,
                 list = TRUE, extra = list(tree,traits))

  return (samples)
}

compute.postquant <- function(samples,initial){

  num.param = ncol(samples)
  num.samples = nrow(samples)

  res = sapply(1:num.param, function(i) (sum(samples[,i]<initial[i])) / num.samples)
  return(res)
}
