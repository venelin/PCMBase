source('postQuantileTools.R')
library(adaptMCMC)
library(PCMBase)
library(mvtnorm)

#1. first read from the command line arguments the replication id (an integer number)

args<-commandArgs(TRUE)

replication = 0
# test if there is exactly one argument: if not, return an error

if (length(args)!=1) {
  stop("One argument must be supplied! It is the replication id", call.=FALSE)

} else {
  # default output file
  replication = as.numeric(args[1])
}

# 2. use this integer to set the random generator seed (set.seed())

set.seed(replication)

# 3. read the tree into a global variable

tree = get(load("Tree.RData"))

# 4. param.replicates <- generate.seeds(replications) (just one replication!)

param.replicates <- generate.seeds(1)

# 5. go on with the replication

Alpha1 = param.replicates[['Alpha1']][[1]]
Alpha2 = param.replicates[['Alpha2']][[1]]
Theta = param.replicates[['Theta']][[1]]
Sigma = param.replicates[['Sigma']][[1]]
Sigmae = param.replicates[['Sigmae']][[1]]

param.m = list('Alpha1' = Alpha1,'Alpha2' = Alpha2,'Theta'= Theta,'Sigma' = Sigma,'Sigmae'=Sigmae)  ######### or I take out the names

model = generate.model(param.m)

data = generate.data(tree,model)
traits = data$traits

param.l = c(Alpha1[1,1],Alpha1[2,2],Alpha2[1,1],Alpha2[2,2],Theta,
            Sigma[1,1],Sigma[1,2],Sigma[2,2],Sigmae[1,1],Sigmae[2,2])

samples_res = generate.samples(posterior,300000,tree,traits)

samp.coda <- convert.to.coda(samples_res)

post.quant = compute.postquant(samples_res$samples,param.l)


# ... after executing the MCMC

# 6. to store the resulting posterior quantile vector in a file job_replication#id.RData

save(post.quant, file=paste("job_replication",replication,".RData",sep=""))
save(samp.coda, file=paste("coda_sample_job_replication",replication,".RData",sep=""))
# 7. once this is tested locally, the next task is to copy everything to the cluster euler:
# 7.1. copy and install the package PCMBase on the cluster
