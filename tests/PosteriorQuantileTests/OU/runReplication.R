source('postQuantileTools.R')
library(adaptMCMC)
library(PCMBase)
library(mvtnorm)
# First read from the command line arguments the replication id (an integer number)

args<-commandArgs(TRUE)

replication = 0
# Test if there is exactly one argument: if not, return an error

if (length(args)!=1) {
  stop("One argument must be supplied! It is the replication id", call.=FALSE)

} else {
  # default output file
  replication = as.numeric(args[1])
}

# Use this integer to set the random generator seed (set.seed())

set.seed(replication)

# Read the tree into a global variable

tree = get(load("Tree.RData"))

# Generate the seeds.

param.replicates <- generate.seeds(1)

# Parse the parameters

Alpha = param.replicates[['Alpha']][[1]]
Theta = param.replicates[['Theta']][[1]]
Sigma = param.replicates[['Sigma']][[1]]
Sigmae = param.replicates[['Sigmae']][[1]]

# Create a list of matrices that includes all parameters

param.m = list('Alpha' = Alpha,'Theta'= Theta,'Sigma' = Sigma,'Sigmae'=Sigmae)  ######### or I take out the names

# Generate the model

model = generate.model(param.m)

# Generate the traits

data = generate.data(tree,model)
traits = data$traits

#Create a list where each entry is the value of a parameter. Required step for MCMC

param.l = c(Alpha[1,1],Alpha[2,2],Theta,Sigma[1,1],Sigma[1,2],Sigma[2,2],Sigmae[1,1],Sigmae[2,2])

# Generate the samples

samples_res = generate.samples(posterior,100,tree,traits)

# Convert the samples to coda objects

samp.coda <- convert.to.coda(samples_res)

# Compute the posterior quantiles for the parameters using the samples generated

post.quant = compute.postquant(samples_res$samples,param.l)

# Store the resulting posterior quantile vector and the coda object

save(post.quant, file=paste("job_replication",replication,".RData",sep=""))
save(samp.coda, file=paste("coda_sample_job_replication",replication,".RData",sep=""))
